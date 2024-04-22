#!/usr/bin/python
"""
This script generates Overhangs by stringing together Arcs, allowing successful fdm-3d-printing of large 90 deg overhangs!
The genius idea emerged from Steven McCulloch, who coded a demonstration and the basic mechanics: https://github.com/stmcculloch/arc-overhang
This python script builds up on that and offers a convenient way to integrate the ArcOverhangs into an existing gcode-file.
HOW TO USE:
Option A) open your system console and type 'python ' followed by the path to this script and the path of the gcode file. This will overwrite the file.
Option B) open PrusaSlicer, go to print-settings-tab -> output-options. Locate the window for post-processing-script.
    In that window enter: the full path to your python exe, an empty space, and the full path to this script.
    If the either path contains any empty spaces, mask them as described here: https://manual.slic3r.org/advanced/post-processing
=>PrusaSlicer will execute the script after the export of the Gcode, therefore the view in the window wont change. Open the finished gcode file to see the results.
If you want to change generation settings: Scroll to 'Parameter' section. Settings from PrusaSlicer will be extracted automatically from the gcode.
Requirements:
Python 3.5+ and the librarys: shapely 1.8+, numpy 1.2+, numpy-hilbert-curve matplotlib for debugging
Slicing in PrusaSlicer is mandatory.
Tested only in PrusaSlicer 2.5&Python 3.10, other versions might need adapted keywords.
Notes:
This code is a little messy. Usually I would divide it into multiple files, but that would compromise the ease of use.
Therefore, I divided the code into sections, marked with ###
Feel free to give it some refactoring and add more functionalities!
Used Coding-Flavour: variable Names: smallStartEveryWordCapitalized, 'to' replaced by '2', same for "for"->"4". Parameters: BigStartEveryWordCapitalized
Known issues:
-pointsPerCircle>80 might give weird results
-MaxDistanceFromPerimeter >=2*perimeterwidth might give a weird result.
-avoid using the code multiple times onto the same gcode, since the bridge infill is deleted when the arcs are generated.
"""
import sys
import json
import time
import copy
from shapely import Point, Polygon, LineString
from shapely.ops import nearest_points
import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval
import warnings
import random
import platform
import argparse
from arc_library import Layer, plot_geometry, Arc, count_lines_in_main
from config import makeFullSettingDict

def plot_arcs_and_geometry(idx, final_arcs, arcs, remaining_space, start_line_string, start_pt):
    """
    Plot arcs and other geometric data for visualization.

    Args:
    idx (int): Current iteration index.
    final_arcs (list): List of final arc geometries.
    arcs (list): List of current arc geometries.
    remaining_space (Geometry): The remaining geometric space.
    start_line_string (Geometry): The starting line string geometry.
    start_pt (Geometry): The starting point geometry.
    plot_params (dict): Dictionary of plotting parameters.
    """
    plt.title(f"Iteration {idx}, Total No Start Points: {len(final_arcs)}, Total No Arcs: {len(arcs)}")
    plot_geometry(start_line_string, 'r')
    plot_geometry([arc.poly for arc in arcs], changecolor=True)
    plot_geometry(remaining_space, 'g', filled=True)
    plot_geometry(start_pt, "r")
    plt.axis('square')
    plt.show()


def getFileStreamAndPath(args, read=True):
    try:
        if read:
            f = open(args.path, "r")
        else:
            f=open(args.path, "w")
        return f,args.path
    except IOError:
        input("File not found.Press enter.")
        sys.exit(1)

def splitGCodeIntoLayers(gcode:list)->list:
    gcode_list = []
    buff=[]
    for linenumber,line in enumerate(gcode):
        if ";LAYER_CHANGE" in line:
            gcode_list.append(buff)
            buff=[]
            buff.append(line)
        else:
            buff.append(line)
    gcode_list.append(buff)  #catch last layer
    print("last read linenumber:",linenumber)
    return gcode_list



################################# HELPER FUNCITONS Polygon->Arc #################################
#################################################################################################

def midpoint(p1:Point, p2:Point):
    return Point((p1.x + p2.x)/2, (p1.y + p2.y)/2)

def getStartPtOnLS(ls:LineString,kwargs:dict={},choseRandom:bool=False)->Point:
    if ls.geom_type=="MultiLineString" or ls.geom_type=="GeometryCollection":
        lengths=[]
        for lss in ls.geoms:
            if lss.geom_type=="LineString":
                lengths.append(lss.length)
            else:
                print("Startline Item bizzare Type of geometry:",lss.geom_type)
                lengths.append(0)
        lsidx=np.argmax(lengths)
        if not lsidx.is_integer():
            try:
                lsidx=lsidx[0]#if multiple max values: take first occurence
            except:
                print("Exception used for lsidx, should be int:", lsidx)
                lsidx=0
        ls=ls.geoms[lsidx]
    if len(ls.coords)<2:
        warnings.warn("Start LineString with <2 Points invalid")
        input("Can not run script, gcode unmodified. Press Enter")
        raise ValueError("Start LineString with <2 Points invalid")
    if len(ls.coords)==2:
        return midpoint(Point(ls.coords[0]),Point(ls.coords[1]))
    scores=[]
    curLength=0
    pts=[Point(p) for p in ls.coords]
    if choseRandom:
        return random.choice(pts)
    coords=[np.array(p) for p in ls.coords]
    for idp,p in enumerate(pts):
        if idp==0 or idp==len(pts)-1:
            scores.append(0)
            continue
        curLength+=p.distance(pts[idp-1])
        relLength=curLength/ls.length
        lengthscore=1-np.abs(relLength-0.5)#Hat-function: pointscore=1 at relLength=0.5, 0 at start or end.
        v1=coords[idp]-coords[idp-1]
        v2=coords[idp+1]-coords[idp]
        if np.linalg.norm(v1)>0 and np.linalg.norm(v2)>0:#calc angle only for non-zero-vectors
            v1=v1/np.linalg.norm(v1)
            v2=v2/np.linalg.norm(v2)
            anglescore=np.abs(np.sin(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))))#prefer points at corners
            anglescore*=kwargs.get("CornerImportanceMultiplier",1)
            scores.append(lengthscore+anglescore)
        else:
            scores.append(lengthscore)
    maxIndex=scores.index(max(scores))
    return pts[maxIndex]


def get_farthest_point(arc, base_poly, remaining_empty_space):
    longest_distance = -1
    farthest_point = Point(0, 0)
    pointFound = False

    # Calculate buffer outside the loop
    buffer_poly = remaining_empty_space.buffer(1e-2)

    if arc.geom_type == 'Polygon':
        arc_coords = np.array(arc.exterior.coords)
    elif arc.geom_type == 'LineString':
        arc_coords = np.linspace(arc.coords[0], arc.coords[1])
    else:
        print('get_farthest_distance: Wrong shape type given', type(arc))
        # Add error handling if needed

    # Bounding box of the base polygon
    min_x, min_y, max_x, max_y = base_poly.bounds

    # Bounding box check
    arc_coords = arc_coords[(arc_coords[:,0] >= min_x) & (arc_coords[:,0] <= max_x) & 
                            (arc_coords[:,1] >= min_y) & (arc_coords[:,1] <= max_y)]

    for p in arc_coords:
        point = Point(p)
        distance = point.distance(base_poly.boundary)
        if distance > longest_distance and buffer_poly.contains(point):
            longest_distance = distance
            farthest_point = point
            pointFound = True

    if pointFound:
        # Use nearest_points from shapely.ops
        point_on_poly = nearest_points(base_poly, farthest_point)[0]
        return farthest_point, longest_distance, point_on_poly
    else:
        return None, None, None

def move_toward_point(start_point:Point, target_point:Point, distance:float)->Point:
    """Moves a point a set distance toward another point"""

    # Calculate the direction in which to move
    dx = target_point.x - start_point.x
    dy = target_point.y - start_point.y

    # Normalize the direction
    magnitude = (dx**2 + dy**2)**0.5
    dx /= magnitude
    dy /= magnitude

    # Move the point in the direction of the target by the set distance
    return Point(start_point.x + dx*distance, start_point.y + dy*distance)

def redistribute_vertices(geom:LineString, distance:float)->LineString:
    if geom.geom_type == 'LineString':
        num_vert = int(round(geom.length / distance))
        if num_vert == 0:
            num_vert = 1
        return LineString(
            [geom.interpolate(float(n) / num_vert, normalized=True)
             for n in range(num_vert + 1)])
    elif geom.geom_type == 'MultiLineString':
        parts = [redistribute_vertices(part, distance) for part in geom.geoms]
        return type(geom)([p for p in parts if not p.is_empty])
    else:
        warnings.warn('unhandled geometry %s', (geom.geom_type,))
        return geom

def generateMultipleConcentricArcs(startpt:Point,rMin:float,rMax:float, boundaryLineString:LineString,remainingSpace:Polygon,kwargs={})->list:
    arcs=[]
    r=rMin
    while r<=rMax:
        arcObj=Arc(startpt,r,kwargs=kwargs)
        arc=arcObj.generateConcentricArc(startpt,remainingSpace)
        if arc.intersects(boundaryLineString) and not kwargs.get("UseLeastAmountOfCenterPoints",False):
            break
        arcs.append(arcObj)
        #print("True Arc type:",type(arc4gcode))
        r+=kwargs.get("ArcWidth")
    return arcs

################################# HELPER FUNCTIONS Arc->GCode #################################
###############################################################################################

def getArcBoundarys(concentricArcs:list)->list:
    '''Handle arcs composited from multiple parts'''
    boundarys=[]
    for arc in concentricArcs:
        arcLine=arc.extractArcBoundary()
        if type(arcLine)==type([]):
            for arc in arcLine:
                boundarys.append(arc.arcline)
        else:
            boundarys.append(arcLine)
    return boundarys

def readSettingsFromGCode2dict(gcodeLines:list,fallbackValuesDict:dict)->dict:
    gCodeSettingDict=fallbackValuesDict
    isSetting=False
    for line in gcodeLines:
        if "; prusaslicer_config = begin" in line or "CONFIG_BLOCK_START" in line:
            isSetting=True
            continue
        if isSetting :
            setting=line.strip(";").strip("\n").split("= ")
            if len(setting)==2:
                try:
                    gCodeSettingDict[setting[0].strip(" ")]=literal_eval(setting[1]) # automaticly convert into int,float,...
                except:
                    gCodeSettingDict[setting[0].strip(" ")]=setting[1] # leave the complex settings as strings. They shall be handled individually if necessary
            elif len(setting)>2:
                gCodeSettingDict[setting[0].strip(" ")]=setting[1:]
                warnings.warn(f"PrusaSlicer Setting {setting[0]} not in the expected key/value format, but added into the settings-dictionarry")
            else:
                print("Could not read setting from PrusaSlicer:",setting)
    if "%" in str(gCodeSettingDict.get("perimeter_extrusion_width")) : #overwrite Percentage width as suggested by 5axes via github
        gCodeSettingDict["perimeter_extrusion_width"]=gCodeSettingDict.get("nozzle_diameter")*(float(gCodeSettingDict.get("perimeter_extrusion_width").strip("%"))/100)
    isWarned=False
    for key,val in gCodeSettingDict.items():
        if isinstance(val,tuple) :
            if gCodeSettingDict.get("Fallback_"+key):
                gCodeSettingDict[key]=gCodeSettingDict.get("Fallback_"+key)
                #warnings.warn(f"{key}: Fallback value used: {gCodeSettingDict.get(key)}")
            else:
                gCodeSettingDict[key]=val[0]
                if not isWarned:
                    if key != 'wiping_volumes_extruders':
                        warnings.warn(f"{key = } {type(val)} was specified as tuple/list, this is normal for using multiple extruders. For all list values First values will be used. If unhappy: Add manual fallback value by searching for ADD FALLBACK in the code. And add 'Fallback_<key>:<yourValue>' into the dictionary.")
                        isWarned=True
    return gCodeSettingDict

def checkforNecesarrySettings(gCodeSettingDict:dict)->bool:
    if not gCodeSettingDict.get("use_relative_e_distances"):
        warnings.warn("Script only works with relative e-distances enabled in PrusaSlicer. Change acordingly.")
        return False
    if gCodeSettingDict.get("extrusion_width")<0.001 or gCodeSettingDict.get("perimeter_extrusion_width")<0.001 or gCodeSettingDict.get("solid_infill_extrusion_width")<0.001:
        warnings.warn("Script only works with extrusion_width and perimeter_extrusion_width and solid_infill_extrusion_width>0. Change in PrusaSlicer acordingly.")
        return False
    if not gCodeSettingDict.get("overhangs"):
        warnings.warn("Overhang detection disabled in PrusaSlicer. Activate in PrusaSlicer for script success!")
        return False
    if gCodeSettingDict.get("bridge_speed")>5:
        warnings.warn(f"Your Bridging Speed is set to {gCodeSettingDict.get('bridge_speed'):.0f} mm/s in PrusaSlicer. This can cause problems with warping.<=5mm/s is recommended")
    if gCodeSettingDict.get("infill_first"):
        warnings.warn("Infill set in PrusaSlicer to be printed before perimeter. This can cause problems with the script.")
    if gCodeSettingDict.get("external_perimeters_first"):
        warnings.warn("PrusaSlicer-Setting: External perimeter is printed before inner perimeters. Change for better overhang performance. ")
    if not gCodeSettingDict.get("avoid_crossing_perimeters"):
        warnings.warn("PrusaSlicer-Setting: Travel Moves may cross the outline and therefore cause artefacts in arc generation.")
    return True

def calcEStepsPerMM(settingsdict:dict,layerheight:float=None)->float:
    if layerheight:# case: printing on surface.
        w=settingsdict.get("infill_extrusion_width")
        h=layerheight
        eVol=(w-h)*h+np.pi*(h/2)**2 *settingsdict.get("HilbertInfillExtrusionMultiplier",1)
    else:   #case: bridging, used for arcs.
        eVol = (settingsdict.get("nozzle_diameter")/2)**2 * np.pi *settingsdict.get("ArcExtrusionMultiplier",1)#printing in midair will result in circular shape. Scource: https://manual.slic3r.org/advanced/flow-math
    if settingsdict.get("use_volumetric_e"):
        return eVol
    else:
        eInMm=eVol/((settingsdict.get("filament_diameter")/2)**2 *np.pi)
        return eInMm

def p2GCode(p:Point,E=0,**kwargs)->str:
    line=f"G1 X{p.x:.6} Y{p.y:.6} "
    line+="E0" if E==0 else f"E{E:.7f}"
    if kwargs.get('F'):
        line+=f" F{kwargs.get('F'):0d}"
    line+='\n'
    return line

def retractGCode(retract:bool=True,kwargs:dict={})->str:
    retractDist=kwargs.get("retract_length",1)
    restartExtra = kwargs.get("retract_restart_extra")
    E= -retractDist if retract else retractDist + restartExtra
    return f"G1 E{E} F{kwargs.get('retract_speed',35)*60}\n"

def setFeedRateGCode(F:int)->str:
    return f"G1 F{F}\n"

def arc2GCode(arcline:LineString,eStepsPerMM:float,arcidx=None,kwargs={})->list:
    GCodeLines=[]
    p1=None
    pts=[Point(p) for p in arcline.coords]
    if len(pts)<2:
        return []
    #plt.plot([p.x for p in pts],[p.y for p in pts])
    #plt.axis('square')
    #plt.show()
    extDist=kwargs.get("ExtendArcDist",0.5)
    pExtend=move_toward_point(pts[-2],pts[-1],extDist)
    # arcPrintSpeed=np.clip(arcline.length/(kwargs.get("ArcSlowDownBelowThisDuration",3))*60,
    #                         kwargs.get("ArcMinPrintSpeed",1*60),kwargs.get('ArcPrintSpeed',2*60)) # *60 bc unit conversion:mm/s=>mm/min
    arcPrintSpeed = kwargs.get('ArcPrintSpeed', 90)
    # print(f'{arcPrintSpeed = }')
    for idp,p in enumerate(pts):
        if idp==0:
            GCodeLines.append(retractGCode(retract=True,kwargs=kwargs))
            p1=p
            GCodeLines.append(f";Arc {arcidx if arcidx else ' '} Length:{arcline.length}\n")
            GCodeLines.append(p2GCode(p,F=kwargs.get('ArcTravelFeedRate',100*60)))#feedrate is mm/min...
            GCodeLines.append(retractGCode(retract=False,kwargs=kwargs))
            GCodeLines.append(setFeedRateGCode(arcPrintSpeed))
        else:
            dist=p.distance(p1)
            if dist>kwargs.get("GCodeArcPtMinDist",0.1):
                GCodeLines.append(p2GCode(p,E=dist*eStepsPerMM))
                p1=p
        if idp==len(pts)-1:
            GCodeLines.append(p2GCode(pExtend,E=extDist*eStepsPerMM))#extend arc tangentially for better bonding between arcs
    return GCodeLines        

# def hilbert2GCode(allhilbertpts:list,parameters:dict,layerheight:float):
#     hilbertGCode=[]
#     eStepsPerMM=calcEStepsPerMM(parameters,layerheight)
#     for idc,curvepts in enumerate(allhilbertpts):
#         for idp,p in enumerate(curvepts):
#             if idp==0:
#                 hilbertGCode.append(p2GCode(p,F=parameters.get("ArcTravelFeedRate")))
#                 if idc==0:
#                     hilbertGCode.append(retractGCode(False,parameters))
#             elif idp==1:
#                 hilbertGCode.append(p2GCode(p,E=eStepsPerMM*p.distance(lastP), F=parameters.get("aboveArcsInfillPrintSpeed")))
#             else:
#                 hilbertGCode.append(p2GCode(p,E=eStepsPerMM*p.distance(lastP)))
#             lastP=p
#         #finish line
#     hilbertGCode.append(retractGCode(True,parameters))
#     return hilbertGCode

def _warning(message,category = UserWarning, filename = '', lineno = -1,*args, **kwargs):
    print(f"{filename}:{lineno}: {message}")
warnings.showwarning = _warning

def main(gCodeFileStream,path2GCode,skipInput,overrideSettings)->None:
    '''Here all the work is done, therefore it is much to long.'''
    gCodeLines=gCodeFileStream.readlines()
    defaultSettings = {
    "Fallback_nozzle_diameter": 0.4,
    "Fallback_filament_diameter": 1.75,
    "Fallback_perimeter_extrusion_width": 0.4} # add fallback settings here

    gCodeSettingDict = readSettingsFromGCode2dict(gCodeLines, defaultSettings)
    parameters = makeFullSettingDict(gCodeSettingDict)
    output_name = path2GCode.replace('.gcode', '-arcs-added.gcode')
    print(f'{output_name = }')
    parameters["Path2Output"] = output_name
    
    if overrideSettings:
        parameters.update(json.loads(overrideSettings))
    if not checkforNecesarrySettings(gCodeSettingDict):
        warnings.warn("Incompatible PursaSlicer-Settings used!")
        input("Can not run script, gcode unmodified. Press enter to close.")
        raise ValueError("Incompatible Settings used!")
    layerobjs=[]; gcodeWasModified=False; numOverhangs=0; lastfansetting=127;
    if not gCodeFileStream:
        print('No file found')
        sys.exit()
    layers=splitGCodeIntoLayers(gCodeLines)
    gCodeFileStream.close()
    print("layers:",len(layers))
    # setup layer objects for layers
    for idl,layerlines in enumerate(layers):
        layer=Layer(layerlines,parameters,idl)
        lastfansetting=layer.spotFanSetting(lastfansetting)
        layerobjs.append(layer)
        
    outputLayers = copy.deepcopy(layerobjs)
    # Iterate through each layer object with its index
    for idl, layer in enumerate(layerobjs):

        # Skip to next iteration if no valid polygons resulted from processing
        # validpolys are overhangs that will be converted to arcs
        # old geom. will be deleted of the poly
        if not layer.validpolys:
            continue
        numOverhangs += 1
        print(f"overhang found layer {idl}:",len(layer.polys), f"Z: {layer.z:.2f}")
        # Prepare to apply special cooling settings based on layer height
        maxZ = layer.z + parameters.get("specialCoolingZdist")
        idoffset = 1 ; arcOverhangGCode = [] 
        currZ = layer.z
    
        # Handling overhangs by printing details and adjusting cooling on subsequent layers
        # Loop through the following layers to apply special settings until the maxZ is reached
        while currZ <= maxZ and idl + idoffset <= len(layerobjs) - 1:
            # Update the current layer's Z height
            currZ = layerobjs[idl + idoffset].z
            # Extend old polygons of subsequent layers with valid polygons from the current layer
            layerobjs[idl + idoffset].oldpolys.extend(layer.validpolys)
            # Move to the next layer
            idoffset += 1
        # Initialize variables for tracking and processing overhangs with arc generation
        prevLayer = layerobjs[idl-1]  # Get the previous layer to work with its geometries
        prevLayer.makeExternalPerimeter2Polys()  # Convert external perimeters to polygons
        
        # Iterate over valid polygons in the current layer to generate arcs
        for poly in layer.validpolys:
            # Retrieve parameters for arc generation
            MaxDistanceFromPerimeter = parameters.get("MaxDistanceFromPerimeter"); rMax = parameters.get("RMax", 15); 
            arcWidth = parameters.get("ArcWidth"); rMinStart = parameters.get("RMinStartMultiple") * parameters.get("nozzle_diameter"); 
            rMin = rMinStart; finalarcs = []; arcs = []; arcs4gcode = []
        
            # Determine starting point for arcs based on the previous layer's perimeter
            startLineString, boundaryWithOutStartLine = prevLayer.makeStartLineString(poly, parameters)
            if startLineString is None:
                warnings.warn("Skipping Polygon because no StartLine Found")
                continue
            startpt = getStartPtOnLS(startLineString, parameters)
            remainingSpace = poly
        
            # Generate concentric arcs from the start point within defined boundaries
            concentricArcs = generateMultipleConcentricArcs(startpt, rMinStart, rMax, boundaryWithOutStartLine, remainingSpace, parameters)
            # Optionally, print the number of arcs generated for debugging
            # print(f"number of concentric arcs generated:", len(concentricArcs))

            # Error handling for inadequate arc generation from the chosen starting point
            attempt_count = 0; max_attempts = 20  # Define a limit to avoid infinite loops
            
            while len(concentricArcs) < parameters.get("MinStartArcs") and attempt_count < max_attempts:
                if attempt_count == 0:
                    # Initial attempt to adjust the start point by redistributing vertices
                    startpt = getStartPtOnLS(redistribute_vertices(startLineString, 0.1), parameters)
                elif attempt_count < 10:
                    # Next 10 attempts: Using random start points
                    print(f"Layer {idl}: Using random Startpoint")
                    startpt = getStartPtOnLS(startLineString, parameters, choseRandom=True)
                else:
                    # Subsequent 10 attempts: Redistribute vertices and choose randomly
                    startpt = getStartPtOnLS(redistribute_vertices(startLineString, 0.1), parameters, choseRandom=True)
                
                concentricArcs = generateMultipleConcentricArcs(startpt, rMinStart, rMax, boundaryWithOutStartLine, remainingSpace, parameters)
                attempt_count += 1
            
            # Issue a warning if all attempts fail to generate a sufficient number of arcs
            if len(concentricArcs) < parameters.get("MinStartArcs"):
                warnings.warn("Initialization Error: no concentric Arc could be generated at startpoints, moving on")
                continue

            # # Error handling for inadequate arc generation from the chosen starting point
            # if len(concentricArcs) < parameters.get("MinStartArcs"):
            #     # Adjust the start point by redistributing vertices and retry arc generation
            #     startpt = getStartPtOnLS(redistribute_vertices(startLineString, 0.1), parameters)
            #     concentricArcs = generateMultipleConcentricArcs(startpt, rMinStart, rMax, boundaryWithOutStartLine, remainingSpace, parameters)
            
            #     # If the minimum number of arcs is still not met, attempt to find a new starting point randomly
            #     if len(concentricArcs) < parameters.get("MinStartArcs"):
            #         print(f"Layer {idl}: Using random Startpoint")
            #         for idr in range(10):  # Try up to 10 different random start points
            #             startpt = getStartPtOnLS(startLineString, parameters, choseRandom=True)
            #             concentricArcs = generateMultipleConcentricArcs(startpt, rMinStart, rMax, boundaryWithOutStartLine, remainingSpace, parameters)
            #             if len(concentricArcs) >= parameters.get("MinStartArcs"):
            #                 break
            
            #         # If random start points also fail, retry with redistributed vertices and random choices
            #         if len(concentricArcs) < parameters.get("MinStartArcs"):
            #             for idr in range(10):
            #                 startpt = getStartPtOnLS(redistribute_vertices(startLineString, 0.1), parameters, choseRandom=True)
            #                 concentricArcs = generateMultipleConcentricArcs(startpt, rMinStart, rMax, boundaryWithOutStartLine, remainingSpace, parameters)
            #                 if len(concentricArcs) >= parameters.get("MinStartArcs"):
            #                     break
            
            #         # Issue a warning if all attempts fail to generate a sufficient number of arcs
            #         if len(concentricArcs) < parameters.get("MinStartArcs"):
            #             warnings.warn("Initialization Error: no concentric Arc could be generated at startpoints, moving on")
            #             continue


            # Append the last concentric arc to the list of final arcs, usually the smallest or innermost
            finalarcs.append(concentricArcs[-1])
            
            # Process each concentric arc to update the remaining space by subtracting the area covered by the arc
            for arc in concentricArcs:
                # Subtract the buffered polygon of the arc from the remaining space to avoid overlap in future arcs
                remainingSpace = remainingSpace.difference(arc.poly.buffer(0))
                arcs.append(arc)
            
            arcBoundarys = getArcBoundarys(concentricArcs)
            # Collect the boundaries of arcs to prepare for generating G-code instructions
            for arcboundary in arcBoundarys:
                arcs4gcode.append(arcboundary)

            # Start a breadth-first search (BFS) to fill the remaining space using arc information
            idx = 0; safetyBreak = 0; triedFixing = False
            
            while idx < len(finalarcs):
                # Console output management to update status in the same line
                sys.stdout.write("\033[F")  # Move back to the previous line
                sys.stdout.write("\033[K")  # Clear the current line
                print("while executed:", idx, len(finalarcs))  # Display current loop execution status
            
                curArc = finalarcs[idx]
                # Determine the farthest point on the arc from remaining space and check for MultiPolygon geometry
                if curArc.poly.geom_type == "MultiPolygon":
                    farthestPointOnArc, longestDistance, NearestPointOnPoly = get_farthest_point(curArc.poly.geoms[0], poly, remainingSpace)
                else:
                    farthestPointOnArc, longestDistance, NearestPointOnPoly = get_farthest_point(curArc.poly, poly, remainingSpace)
            
                # Check if any viable points are left on the arc; if not, move to the next one
                if not farthestPointOnArc or longestDistance < MaxDistanceFromPerimeter:
                    idx += 1  # Go to next arc
                    continue
            
                # Adjust the start point towards the arc center to create new concentric arcs
                startpt = move_toward_point(farthestPointOnArc, curArc.center, parameters.get("ArcCenterOffset", 2))
                concentricArcs = generateMultipleConcentricArcs(startpt, rMin, rMax, poly.boundary, remainingSpace, parameters)
                arcBoundarys = getArcBoundarys(concentricArcs)
            
                # If concentric arcs are generated successfully, update structures and remaining space
                if len(concentricArcs) > 0:
                    for arc in concentricArcs:
                        remainingSpace = remainingSpace.difference(arc.poly.buffer(1e-2))
                        arcs.append(arc)
                    finalarcs.append(concentricArcs[-1])
                    for arcboundary in arcBoundarys:
                        arcs4gcode.append(arcboundary)
                else:
                    idx += 1  # No possible concentric arcs found; arc is complete, proceed to next
            
                # Implement a safety break to avoid infinite loops
                safetyBreak += 1
                if safetyBreak > parameters.get("SafetyBreak_MaxArcNumber", 2000):
                    print('\nHit safety break.')
                    break

                if parameters.get("plotArcsEachStep"):
                    plot_arcs_and_geometry(idx, finalarcs, arcs, remainingSpace, startLineString, startpt)
                
                # Check conditions to identify if arc generation is stuck at the beginning
                if len(finalarcs) == 1 and idx == 1 and remainingSpace.area / poly.area * 100 > 50 and not triedFixing:
                    # Automated fix for situations where arc generation is stuck at a tight spot during the initial steps
                    # Reset parameters and index for arc generation retry
                    parameters["ArcCenterOffset"] = 0; rMin = arcWidth / 1.5; idx = 0; triedFixing = True  
                    print("the arc-generation got stuck at a tight spot during startup. Used Automated fix: set ArcCenterOffset to 0")
                
                # If a fix was attempted and still only one arc has been generated and the problem persists, notify failure
                if triedFixing and len(finalarcs) == 1 and idx == 1:
                    print("fix did not work.")  # Print failure message if the automated fix did not resolve the issue
                
                # Indicate that processing for the current polygon is complete


            remain2FillPercent=remainingSpace.area/poly.area*100
            if  remain2FillPercent> 100-parameters.get("WarnBelowThisFillingPercentage"):
                warnings.warn(f"layer {idl}: The Overhang Area is only {100-remain2FillPercent:.0f}% filled with Arcs. Please try again with adapted Parameters: set 'ExtendIntoPerimeter' higher to enlargen small areas. lower the MaxDistanceFromPerimeter to follow the curvature more precise. Set 'ArcCenterOffset' to 0 to reach delicate areas. ")

            if parameters.get("plotArcsFinal"):
                plot_arcs_and_geometry(idx, finalarcs, arcs, remainingSpace, startLineString, startpt)

            # Calculate extruder steps per millimeter based on given parameters
            eStepsPerMM = calcEStepsPerMM(parameters)
            # Add a command to turn on the cooling fan to the specified speed for bridge settings
            arcOverhangGCode.append(f"M106 S{round(parameters.get('bridge_fan_speed',100)*2.55)}\n")
            
            # Filter out any empty geometries before generating G-code
            arcs4gcode = [x for x in arcs4gcode if not x.is_empty]
            
            # Iterate over the arcs to generate G-code for each
            for ida, arc in enumerate(arcs4gcode):
                # Generate G-code for each arc and append to the list
                arcGCode = arc2GCode(arcline=arc, eStepsPerMM=eStepsPerMM, arcidx=ida, kwargs=parameters)
                arcOverhangGCode.append(arcGCode)
                
                # Optionally add a command for time-lapse photography after every specified number of arcs
                if parameters.get("TimeLapseEveryNArcs") > 0:
                    if ida % parameters.get("TimeLapseEveryNArcs") == 0:
                        arcOverhangGCode.append("M240\n")  # Command to take a photo
            
            # Indicate that modifications have been made to the G-code and that the G-code was successfully modified
            modify = True; gcodeWasModified = True
                    
                            # #apply special cooling settings:
                            # if len(layer.oldpolys)>0 and gcodeWasModified:
                            #     modify=True
                            #     print("oldpolys found in layer:",idl)
                            #     layer.spotSolidInfill()
                            #     layer.makePolysFromSolidInfill(extend=parameters.get("ExtendIntoPerimeter"))
                            #     layer.solidPolys=layer.mergePolys(layer.solidPolys)
                            #     allhilbertpts=[]
                            #     for poly in layer.solidPolys:
                            #         hilbertpts=layer.createHilbertCurveInPoly(poly)
                            #         allhilbertpts.extend(hilbertpts)
                            #         if parameters.get("plotEachHilbert"):
                            #             plot_geometry(hilbertpts,changecolor=True)
                            #             plot_geometry(layer.solidPolys)
                            #             plt.title("Debug")
                            #             plt.axis('square')
                            #             plt.show()
                       # Check if modifications have been flagged for the current layer
        if modify:
            # Create a new Layer instance with empty features and other necessary parameters
            modifiedlayer = Layer([], parameters, idl)  
            # Initialize state variables for G-code injection
            isInjected = False; curPrintSpeed = "G1 F600"; messedWithSpeed = False; messedWithFan = False  

            # If G-code was flagged as modified, prepare to update the layer's G-code
            if gcodeWasModified:
                # Remove the 'Bridge' features from the layer based on valid polygons
                layer.prepareDeletion(featurename="Bridge", polys=layer.validpolys)
                # Check if there are old polygons and prepare to remove 'Solid' features
                if len(layer.oldpolys) > 0:
                    layer.prepareDeletion(featurename=":Solid", polys=layer.oldpolys)
            
            # Placeholder for tracking the start of the G-code injection
            injectionStart = None
            print("modifying GCode")
                    
             # Iterate through each line of the current layer's G-code
            for idline, line in enumerate(layer.lines):
                # Skip processing if there are no valid polygons in this layer
                if not layer.validpolys:
                    continue
                
                # Look for a line containing the type of G-code section, to inject arc G-code at the start of the section
                if ";TYPE" in line and not isInjected:
                    injectionStart = idline  # Mark the injection start point
                    modifiedlayer.lines.append(";TYPE:Arc infill\n")  # Indicate the start of arc infill section
                    modifiedlayer.lines.append(f"M106 S{parameters.get('ArcFanSpeed')}\n")  # Set fan speed for arc printing
                    
                    # Loop through each arc G-code sequence to be injected
                    for overhangline in arcOverhangGCode:
                        for arcline in overhangline:
                            for cmdline in arcline:
                                modifiedlayer.lines.append(cmdline)  # Append each command line to the modified layer's G-code
                    
                    isInjected = True  # Set flag indicating that arc G-code has been injected
            
                    # Add G-code to restore the previous tool position before the injection
                    for line2Check in reversed(layer.lines[:injectionStart]):
                        if "X" in line2Check:
                            modifiedlayer.lines.append(line2Check)  
                            break

                                                        # if layer.oldpolys:
                                                        #     if ";TYPE" in line and not hilbertIsInjected:# startpoint of solid infill: print all hilberts from here.
                                                        #         hilbertIsInjected=True
                                                        #         injectionStart=idline
                                                        #         modifiedlayer.lines.append(";TYPE:Solid infill\n")
                                                        #         modifiedlayer.lines.append(f"M106 S{parameters.get('aboveArcsFanSpeed')}\n")
                                                        #         hilbertGCode=hilbert2GCode(allhilbertpts,parameters,layer.height)
                                                        #         modifiedlayer.lines.extend(hilbertGCode)
                                                        #         #add restored pre-injected tool position
                                                        #         for id in reversed(range(injectionStart)):
                                                        #             if "X" in layer.lines[id]:
                                                        #                 modifiedlayer.lines.append(layer.lines[id])
                                                        #                 break
                        # Check for specific G-code commands to adjust settings based on line content
                        if "G1 F" in line.split(";")[0]:  # Check if the line contains a speed command, ignoring comments
                            curPrintSpeed = line  # Update the current print speed setting
                
                # Evaluate whether the current line should be exported based on layer conditions
                if layer.exportThisLine(idline):
                    # Adjust print settings if the line is close to bridging areas
                    if layer.isClose2Bridging(line, parameters.get("CoolingSettingDetectionDistance")):
                        # Modify fan speed for optimal cooling during bridging, if not already done
                        if not messedWithFan:
                            modifiedlayer.lines.append(f"M106 S{parameters.get('aboveArcsFanSpeed')}\n")
                            messedWithFan = True
                
                        # Append modified line with updated print speed for areas above arcs
                        modline = line.strip("\n") + f" F{parameters.get('aboveArcsPerimeterPrintSpeed')}\n"
                        modifiedlayer.lines.append(modline)
                        messedWithSpeed = True
                    else:
                        # Reset fan speed to original if specific layer-wide settings are not applied
                        if messedWithFan and not parameters.get("applyAboveFanSpeedToWholeLayer"):
                            modifiedlayer.lines.append(f"M106 S{layer.fansetting:.0f}\n")
                            messedWithFan = False
                        
                        # Revert to the original print speed if it was changed
                        if messedWithSpeed:
                            modifiedlayer.lines.append(curPrintSpeed + "\n")
                            messedWithSpeed = False
                        
                        # Append the original line to the modified layer after any necessary adjustments
                        modifiedlayer.lines.append(line)
                        
                        
            if messedWithFan:
                modifiedlayer.lines.append(f"M106 S{layer.fansetting:.0f}\n")
                messedWithFan=False
            outputLayers[idl]=modifiedlayer  # overwrite the infos

    if gcodeWasModified:
        path2GCode=parameters.get("Path2Output")
        f=open(path2GCode,"w")
        print("write to",path2GCode)
        for layer in outputLayers:
            f.writelines(layer.lines)
        print('\nAdding settings')
        f.write('\n---------Arc Overhang Settings-----------')
        for k,v in parameters.items():
            f.write(f'\n;{k}: {v}')
        f.close()
    
    else:
        if numOverhangs > 0:
            print(f"Found {numOverhangs} overhangs, but no arcs could be generated due to unusual geometry.")
        else:
            print(f"Analysed {len(layerobjs)} Layers, but no matching overhangs found->no arcs generated. If unexpected: look if restricting settings like 'minArea' or 'MinBridgeLength' are correct.")
    print("Script execution complete.")
    if not skipInput:
        input("Press enter to exit.")

def parse_args():
    parser = argparse.ArgumentParser(description="Process G-code files.")
    parser.add_argument('path', type=str, help='Path to the G-code file')
    parser.add_argument('--settings', type=str, help='JSON string of settings to override default values')
    parser.add_argument('--skip-input', action='store_true', help='Skip any user input prompts (non-Windows only)')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    overrideSettings = args.settings if args.settings else None
    # Get file stream and path based on the provided path
    gCodeFileStream, path2GCode = getFileStreamAndPath(args)
    # Determine whether to skip input based on the platform and command line argument
    skipInput = args.skip_input or platform.system() != "Windows"
    
    # Start timing
    start_time = time.perf_counter()
    
    # Call the main function
    main(gCodeFileStream, path2GCode, skipInput, overrideSettings)
    
    # Calculate execution time
    end_time = time.perf_counter()
    execution_time = end_time - start_time
    
    print(f"Execution time: {execution_time} seconds")
    print(f'{count_lines_in_main() = }')
