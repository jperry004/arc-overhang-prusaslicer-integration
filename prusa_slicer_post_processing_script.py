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
#!/usr/bin/python
import sys
import json
from shapely import Point, Polygon, LineString
from shapely.ops import nearest_points
import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval
import warnings
import random
import platform
import argparse
from arc_library import Layer, plot_geometry, Arc
from config import makeFullSettingDict



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

def arc2GCode(arcline:LineString,eStepsPerMM:float,arcidx=None,final_arc=None,kwargs={})->list:
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
    print(f'{arcPrintSpeed = }')
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
    layerobjs=[]
    gcodeWasModified=False
    numOverhangs=0
    if not gCodeFileStream:
        pass
    
    layers=splitGCodeIntoLayers(gCodeLines)
    gCodeFileStream.close()
    print("layers:",len(layers))
    lastfansetting=0 # initialize variable
    for idl,layerlines in enumerate(layers):
        layer=Layer(layerlines,parameters,idl)
        layer.addZ()
        layer.addHeight()
        lastfansetting=layer.spotFanSetting(lastfansetting)
        layerobjs.append(layer)
    for idl,layer in enumerate(layerobjs):
        modify=False
        if idl<1:
            continue # no overhangs in the first layer and dont mess with the setup
    
        layer.extract_features()
        layer.spotBridgeInfill()
        layer.makePolysFromBridgeInfill(extend=parameters.get("ExtendIntoPerimeter",1))
        layer.polys=layer.mergePolys()
        layer.verifyinfillpolys()

        #ARC GENERATION
        if not layer.validpolys:
            continue
        numOverhangs += 1
        print(f"overhang found layer {idl}:",len(layer.polys), f"Z: {layer.z:.2f}")
        #set special cooling settings for the follow up layers
        maxZ=layer.z+parameters.get("specialCoolingZdist")
        idoffset=1
        currZ=layer.z
        while currZ<=maxZ and idl+idoffset<=len(layerobjs)-1:
            currZ=layerobjs[idl+idoffset].z
            layerobjs[idl+idoffset].oldpolys.extend(layer.validpolys)
            idoffset+=1

        #make Startpoint form previous layer
        prevLayer=layerobjs[idl-1]
        prevLayer.makeExternalPerimeter2Polys()
        arcOverhangGCode=[]
        for poly in layer.validpolys:
            #make parameters more readable
            MaxDistanceFromPerimeter=parameters.get("MaxDistanceFromPerimeter") # how much 'bumpiness' you accept in the outline. Lower will generate more small arcs to follow the perimeter better (corners!). Good practice: 2 perimeters+ threshold of 2width=minimal exact touching (if rMin satisfied)
            rMax=parameters.get("RMax",15)
            # pointsPerCircle=parameters.get("PointsPerCircle",80)
            arcWidth=parameters.get("ArcWidth")
            # rMin=parameters.get("ArcCenterOffset")+arcWidth/1.5

            rMinStart = parameters.get("RMinStartMultiple") * parameters.get("nozzle_diameter")
            rMin = rMinStart
            #initialize
            finalarcs=[]
            arcs=[]
            arcs4gcode=[]
            #find StartPoint and StartLineString
            startLineString,boundaryWithOutStartLine=prevLayer.makeStartLineString(poly,parameters)
            if startLineString is None:
                warnings.warn("Skipping Polygon because no StartLine Found")
                continue
            startpt=getStartPtOnLS(startLineString,parameters)
            remainingSpace=poly
            #plot_geometry(thresholdedpoly)
            #plot_geometry(startLineString,'m')
            #plot_geometry(startpt,'r')
            #plt.axis('square')
            #plt.show()
            #first step in Arc Generation

            concentricArcs=generateMultipleConcentricArcs(startpt,rMinStart,rMax,boundaryWithOutStartLine,remainingSpace,parameters)
            #print(f"number of concentric arcs generated:",len(concentricArcs))
            if len(concentricArcs)<parameters.get("MinStartArcs"):
                #possibly bad chosen startpt, errorhandling:
                startpt=getStartPtOnLS(redistribute_vertices(startLineString,0.1),parameters)
                concentricArcs=generateMultipleConcentricArcs(startpt,rMinStart,rMax,boundaryWithOutStartLine,remainingSpace,parameters)
                if len(concentricArcs)<parameters.get("MinStartArcs"):#still insuff start: try random
                    print(f"Layer {idl}: Using random Startpoint")
                    for idr in range(10):
                        startpt=getStartPtOnLS(startLineString,parameters,choseRandom=True)
                        concentricArcs=generateMultipleConcentricArcs(startpt,rMinStart,rMax,boundaryWithOutStartLine,remainingSpace,parameters)
                        if len(concentricArcs)>=parameters.get("MinStartArcs"):
                            break
                    if len(concentricArcs)<parameters.get("MinStartArcs"):
                        for idr in range(10):
                            startpt=getStartPtOnLS(redistribute_vertices(startLineString,0.1),parameters,choseRandom=True)
                            concentricArcs=generateMultipleConcentricArcs(startpt,rMinStart,rMax,boundaryWithOutStartLine,remainingSpace,parameters)
                            if len(concentricArcs)>=parameters.get("MinStartArcs"):
                                break
                    if len(concentricArcs)<parameters.get("MinStartArcs"):
                        warnings.warn("Initialization Error: no concentric Arc could be generated at startpoints, moving on")
                        continue
            arcBoundarys=getArcBoundarys(concentricArcs)
            finalarcs.append(concentricArcs[-1])
            for arc in concentricArcs:
                remainingSpace=remainingSpace.difference(arc.poly.buffer(1e-2))
                arcs.append(arc)
            for arcboundary in arcBoundarys:
                arcs4gcode.append(arcboundary)

            #start bfs (breadth first search algorithm) to fill the remainingspace
            idx=0
            safetyBreak=0
            triedFixing=False
            while idx<len(finalarcs):
                sys.stdout.write("\033[F") #back to previous line
                sys.stdout.write("\033[K") #clear line
                print("while executed:",idx, len(finalarcs))#\r=Cursor at linestart
                curArc=finalarcs[idx]
                if curArc.poly.geom_type=="MultiPolygon":
                    farthestPointOnArc,longestDistance,NearestPointOnPoly=get_farthest_point(curArc.poly.geoms[0],poly,remainingSpace)
                else:
                    farthestPointOnArc,longestDistance,NearestPointOnPoly=get_farthest_point(curArc.poly,poly,remainingSpace)
                if not farthestPointOnArc or longestDistance<MaxDistanceFromPerimeter:#no more pts on arc
                    idx+=1 #go to next arc
                    continue
                startpt=move_toward_point(farthestPointOnArc,curArc.center,parameters.get("ArcCenterOffset",2))
                concentricArcs=generateMultipleConcentricArcs(startpt,rMin,rMax,poly.boundary,remainingSpace,parameters)
                arcBoundarys=getArcBoundarys(concentricArcs)
                #print(f"number of concentric arcs generated:",len(concentricArcs))
                if len(concentricArcs)>0:
                    for arc in concentricArcs:
                        remainingSpace=remainingSpace.difference(arc.poly.buffer(1e-2))
                        arcs.append(arc)
                    finalarcs.append(concentricArcs[-1])
                    for arcboundary in arcBoundarys:
                        arcs4gcode.append(arcboundary)
                else:
                    idx+=1 # no possible concentric arcs found= arc complete, proceed to next
                safetyBreak+=1
                if safetyBreak>parameters.get("SafetyBreak_MaxArcNumber",2000):
                    break
                if parameters.get("plotArcsEachStep"):
                    plt.title(f"Iteration {idx}, Total No Start Points: {len(finalarcs)}, Total No Arcs: {len(arcs)}")
                    plot_geometry(startLineString,'r')
                    plot_geometry([arc.poly for arc in arcs],changecolor=True)
                    plot_geometry(remainingSpace,'g',filled=True)
                    plot_geometry(startpt,"r")
                    plt.axis('square')
                    plt.show()

                if len(finalarcs)==1 and idx==1 and remainingSpace.area/poly.area*100>50 and not triedFixing:
                    #error handling: the arc-generation got stuck at a thight spot during startup. Automated fix:
                    parameters["ArcCenterOffset"]=0
                    rMin=arcWidth/1.5
                    idx=0
                    triedFixing=True
                    print("the arc-generation got stuck at a thight spot during startup. Used Automated fix:set ArcCenterOffset to 0")
                if triedFixing and len(finalarcs)==1 and idx==1:
                    print("fix did not work.")
            #poly finished
            remain2FillPercent=remainingSpace.area/poly.area*100
            if  remain2FillPercent> 100-parameters.get("WarnBelowThisFillingPercentage"):
                warnings.warn(f"layer {idl}: The Overhang Area is only {100-remain2FillPercent:.0f}% filled with Arcs. Please try again with adapted Parameters: set 'ExtendIntoPerimeter' higher to enlargen small areas. lower the MaxDistanceFromPerimeter to follow the curvature more precise. Set 'ArcCenterOffset' to 0 to reach delicate areas. ")
            if parameters.get("plotArcsFinal"):
                plt.title(f"Iteration {idx}, Total No Start Points: {len(finalarcs)}, Total No Arcs: {len(arcs)}")
                plot_geometry(startLineString,'r')
                plot_geometry([arc.poly for arc in arcs],changecolor=True)
                plot_geometry(remainingSpace,'g',filled=True)
                plot_geometry(startpt,"r")
                plt.axis('square')
                plt.show()
            #generate gcode for arc and insert at the beginning of the layer
            eStepsPerMM=calcEStepsPerMM(parameters)
            arcOverhangGCode.append(f"M106 S{np.round(parameters.get('bridge_fan_speed',100)*2.55)}\n")#turn cooling Fan on at Bridge Setting
            #for arc in arcs4gcode:
            #    plot_geometry(arc)
            #    plot_geometry(Point(arc.coords[0]))
            #plt.axis('square')
            #plt.show()
            arcs4gcode = [x for x in arcs4gcode if not x.is_empty]
            for ida,arc in enumerate(arcs4gcode):
                final_arc = ida == len(arcs4gcode) - 1
                if not arc.is_empty:

                    arcGCode=arc2GCode(arcline=arc,eStepsPerMM=eStepsPerMM,arcidx=ida,final_arc=final_arc,kwargs=parameters)
                    arcOverhangGCode.append(arcGCode)
                    if parameters.get("TimeLapseEveryNArcs")>0:
                        if ida%parameters.get("TimeLapseEveryNArcs"):
                            arcOverhangGCode.append("M240\n")

            modify=True
            gcodeWasModified=True

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
        if modify:
            modifiedlayer=Layer([],parameters,idl) # copy the other infos if needed: future to do
            isInjected=False
            # hilbertIsInjected=False
            curPrintSpeed="G1 F600"
            messedWithSpeed=False
            messedWithFan=False
            if gcodeWasModified:
                layer.prepareDeletion(featurename="Bridge",polys=layer.validpolys)
                if len(layer.oldpolys)>0:
                    layer.prepareDeletion(featurename=":Solid",polys=layer.oldpolys)
            #print("FEATURES:",[(f[0],f[2]) for f in layer.features])
            injectionStart=None
            print("modifying GCode")
            for idline,line in enumerate(layer.lines):
                if not layer.validpolys: continue
            
                if ";TYPE" in line and not isInjected:#inject arcs at the very start
                    injectionStart=idline
                    modifiedlayer.lines.append(";TYPE:Arc infill\n")
                    modifiedlayer.lines.append(f"M106 S{parameters.get('ArcFanSpeed')}\n")
                    for overhangline in arcOverhangGCode:
                        for arcline in overhangline:
                            for cmdline in arcline:
                                modifiedlayer.lines.append(cmdline)
                    isInjected=True
                    #add restored pre-injected tool position
                    for id in reversed(range(injectionStart)):
                        if "X" in layer.lines[id]:
                            modifiedlayer.lines.append(layer.lines[id])
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
                if "G1 F" in line.split(";")[0]:#special block-speed-command
                    curPrintSpeed=line
                if layer.exportThisLine(idline):
                    if layer.isClose2Bridging(line,parameters.get("CoolingSettingDetectionDistance")):
                        if not messedWithFan:
                            modifiedlayer.lines.append(f"M106 S{parameters.get('aboveArcsFanSpeed')}\n")
                            messedWithFan=True
                        modline=line.strip("\n")+ f" F{parameters.get('aboveArcsPerimeterPrintSpeed')}\n"
                        modifiedlayer.lines.append(modline)
                        messedWithSpeed=True
                    else:
                        if messedWithFan and not parameters.get("applyAboveFanSpeedToWholeLayer"):
                            modifiedlayer.lines.append(f"M106 S{layer.fansetting:.0f}\n")
                            messedWithFan=False
                        if messedWithSpeed:
                            modifiedlayer.lines.append(curPrintSpeed+"\n")
                            messedWithSpeed=False
                        modifiedlayer.lines.append(line)
            if messedWithFan:
                modifiedlayer.lines.append(f"M106 S{layer.fansetting:.0f}\n")
                messedWithFan=False
            layerobjs[idl]=modifiedlayer  # overwrite the infos

    if gcodeWasModified:
        path2GCode=parameters.get("Path2Output")
        f=open(path2GCode,"w")
        print("write to",path2GCode)
        for layer in layerobjs:
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
    # Call the main function with the arguments
    import time
    
    # Start timing
    start_time = time.perf_counter()
    
    # Call the main function
    main(gCodeFileStream, path2GCode, skipInput, overrideSettings)
    
    # Calculate execution time
    end_time = time.perf_counter()
    execution_time = end_time - start_time
    
    print(f"Execution time: {execution_time} seconds")

