from shapely import Point, Polygon, LineString, GeometryCollection, MultiLineString, MultiPolygon
import warnings
import matplotlib.pyplot as plt
from shapely.ops import linemerge, unary_union
import numpy as np
import random
from hilbert import decode, encode

def count_lines_in_main():
    filename = 'prusa_slicer_post_processing_script.py'
    with open(filename, 'r') as file:
        in_main = False
        main_count = 0

        for line in file:
            stripped_line = line.strip()

            # Check if we are entering the main function
            if stripped_line.startswith("def main("):
                in_main = True
                continue  # Skip the function declaration line

            # Check if we are exiting the main function upon encountering another function definition
            if in_main and stripped_line.startswith("def "):
                break  # Stop counting when another function definition is encountered

            # Count lines only when we are within the main function
            if in_main:
                if stripped_line == "" or stripped_line.startswith("#"):
                    continue  # Skip empty lines and comments within the main function
                main_count += 1

        return main_count




def create_circle(p:Point, radius:float, n:int)->Polygon:
    x=p.x
    y=p.y
    return Polygon([[radius*np.sin(theta)+x, radius*np.cos(theta)+y] for theta in np.linspace(0, 2*np.pi - 2*np.pi/n, int(n))])

def getValueBasedColor(val:float, max_val=10)->tuple:
    normalizedVal = val / max_val
    rgb=[0,0,0]
    rgb[0]=min(normalizedVal,1)
    rgb[2]=1-rgb[0]
    return tuple(rgb)


def plot_geometry(geometry, color='black', linewidth=1,**kwargs):
    if type(geometry)==type([]):
        for idx,geo in enumerate(geometry):
            if kwargs.get("changecolor"):
                color=getValueBasedColor(idx,len(geometry))
            plot_geometry(geo,color=color,linewidth=linewidth,kwargs=kwargs)
    elif geometry.geom_type == 'Point':
        x, y = geometry.x, geometry.y
        plt.scatter(x, y, color=color, linewidth=linewidth)
    elif geometry.geom_type == 'LineString':
        x, y = geometry.xy
        plt.plot(x, y, color=color, linewidth=linewidth)
    elif geometry.geom_type == 'Polygon':
        x, y = geometry.exterior.xy
        plt.plot(x, y, color=color, linewidth=linewidth)
        if kwargs.get("filled"):
            plt.fill(x,y,color=color,alpha=0.8)
        for interior in geometry.interiors:
            x, y = interior.xy
            plt.plot(x, y, color=color, linewidth=linewidth)
            if kwargs.get("filled_holes"):
                plt.fill(x,y,color=color,alpha=0.5)
    elif geometry.geom_type == 'MultiLineString':
        for line in geometry.geoms:
            x, y = line.xy
            plt.plot(x, y, color=color, linewidth=linewidth)
    elif geometry.geom_type == 'MultiPolygon' or geometry.geom_type == "GeometryCollection":
        for polygon in geometry.geoms:
            plot_geometry(polygon,color=color,linewidth=linewidth,kwargs=kwargs)
    else:
        print('Unhandled geometry type: ' + geometry.geom_type)

def p2GCode(p:Point,E=0,**kwargs)->str:
    line=f"G1 X{p.x:.6} Y{p.y:.6} "
    line+="E0" if E==0 else f"E{E:.7f}"
    if kwargs.get('F'):
        line+=f" F{kwargs.get('F'):0d}"
    line+='\n'
    return line

def getPtfromCmd(line:str)->Point:
    x=None
    y=None
    line=line.split(";")[0]
    cmds=line.split(" ")
    for c in cmds:
        if "X" in c:
            x=float(c[1:])
        elif "Y" in c:
            y=float(c[1:])
    if (x is not None) and (y is not None):
        p=Point(x,y)
    else:
        p=None
    return p

def makePolygonFromGCode(lines:list)->Polygon:
    pts=[]
    wiping=False
    for line in lines:
        if line.startswith(";WIPE_END"):
            wiping=False
        elif line.startswith(";WIPE_START"):
            wiping=True

        if wiping:
            continue

        if line.startswith("G1 X"):
            p=getPtfromCmd(line)
            if p:
                pts.append(p)
    if len(pts)>2:
        return Polygon(pts)
    else:
        #print("invalid poly: not enough pts")
        return None

class Layer():
    def __init__(self,lines:list=[],kwargs:dict={},layernumber:int=-1)->None:
        self.lines=lines
        self.layernumber=layernumber
        self.z=kwargs.get("z",None)
        self.polys=[]
        self.validpolys=[]
        self.extPerimeterPolys=[]
        self.binfills=[]
        self.features=[]
        self.oldpolys=[]
        self.dontPerformPerimeterCheck=kwargs.get('notPerformPerimeterCheck',False)
        self.deleteTheseInfills=[]
        self.deletelines=[]
        self.associatedIDs=[]
        self.sinfills=[]
        self.parameters=kwargs
        self.lastP=None
    def extract_features(self)->None:
        buff=[]
        currenttype=""
        start=0
        for idl,line in enumerate(self.lines):
            if ";TYPE:" in line:
                if currenttype:
                    self.features.append([currenttype,buff,start])
                    buff=[]
                    start=idl

                    if self.lines[start - 1].startswith("G1 E") or not "E" in self.lines[start - 1]: # if a travel move came right before this feature:
                        start -= 1
                        while not "E" in self.lines[start - 1]: # include this travel move in the feature (excluding the first retraction step, if applicable)
                            start -= 1
                currenttype=line
            else:
                buff.append(line)
        self.features.append([currenttype,buff,start])# fetch last one
    def addZ(self,z:float=None)->None:
        if z:
            self.z=z
        else:
            for l in self.lines:
                cmd=l.split(";")[0] # work only on the command itself
                if "G1" in cmd and "Z" in cmd:
                    cmds=cmd.split(" ")
                    for c in cmds:
                        if "Z" in c:
                            self.z=float(c[1:])
                            return
    def addHeight(self):
        for l in self.lines:
            if ";HEIGHT" in l:
                h=l.split(":")
                self.height=float(h[-1])
                return
        if self.layernumber > 0:
            warnings.warn(f"Layer {self.layernumber}: no height found, using layerheight default!")
        self.height=self.parameters.get("layer_height")
    def getRealFeatureStartPoint(self,idf:int)->Point:
        """ since GCode only stores destination of the move, the origin of the first move has to be included."""
        if idf<1:
            return None
        lines=self.features[idf-1][1]
        for line in reversed(lines):
            if line.startswith("G1 X"):
                return getPtfromCmd(line)

    def makeExternalPerimeter2Polys(self)->None:
        extPerimeterIsStarted=False
        for idf,fe in enumerate(self.features):
            ftype=fe[0]
            lines=fe[1]

            if "External" in ftype or ("Overhang" in ftype and extPerimeterIsStarted) or ("Overhang" in ftype and self.dontPerformPerimeterCheck): #two different types of perimeter to for a poly: external perimeter and overhang perimeter + option for manual errorhandling, when there is no feature "external"
                if not extPerimeterIsStarted:
                    linesWithStart=[]
                    if idf>1:
                        pt=self.getRealFeatureStartPoint(idf)
                        if type(pt)==type(Point):
                            linesWithStart.append(p2GCode(pt))
                        else:
                            warnings.warn(f"Layer {self.layernumber}: Could not fetch real StartPoint.")
                linesWithStart=linesWithStart+lines
                extPerimeterIsStarted=True
            if extPerimeterIsStarted and (idf==len(self.features)-1 or not ("External" in ftype or "Overhang" in ftype)) :#finish the poly if end of featurelist or different feature
                poly=makePolygonFromGCode(linesWithStart)
                if poly:
                    self.extPerimeterPolys.append(poly)
                extPerimeterIsStarted=False

    def makeStartLineString(self,poly:Polygon,kwargs:dict={}):
        if not self.extPerimeterPolys:
            self.makeExternalPerimeter2Polys()
        if len(self.extPerimeterPolys)<1:
            warnings.warn(f"Layer {self.layernumber}: No ExternalPerimeterPolys found in prev Layer")
            return None,None
        for ep in self.extPerimeterPolys:
            ep=ep.buffer(1e-2)# avoid self intersection error
            if ep.intersects(poly):
                startArea=ep.intersection(poly)
                startLineString=startArea.boundary.intersection(poly.boundary.buffer(1e-2))
                if startLineString.is_empty:
                    if poly.contains(startArea):#if inside no boundarys can overlap.
                        startLineString=startArea.boundary
                        boundaryLineString=poly.boundary
                        if startLineString.is_empty:#still empty? unlikely to happen
                            plt.title("StartLineString is None")
                            plot_geometry(poly,'b')
                            plot_geometry(startArea,filled=True)
                            plot_geometry([ep for ep in self.extPerimeterPolys])
                            plt.legend(["currentLayerPoly","StartArea","prevLayerPoly"])
                            plt.axis('square')
                            plt.show()
                            warnings.warn(f"Layer {self.layernumber}: No Intersection in Boundary,Poly+ExternalPoly")
                            return None,None
                else:
                    boundaryLineString=poly.boundary.difference(startArea.boundary.buffer(1e-2))
                #print("STARTLINESTRING TYPE:",startLineString.geom_type)
                if kwargs.get("plotStart"):
                    print("Geom-Type:",poly.geom_type)
                    plot_geometry(poly,color="b")
                    plot_geometry(ep,'g')
                    plot_geometry(startLineString,color="m")
                    plt.title("Start-Geometry")
                    plt.legend(["Poly4ArcOverhang","External Perimeter prev Layer","StartLine for Arc Generation"])
                    plt.axis('square')
                    plt.show()
                return startLineString,boundaryLineString
        #end of for loop, and no intersection found
        plt.title("no intersection with prev Layer Boundary")
        plot_geometry(poly,'b')
        plot_geometry([ep for ep in self.extPerimeterPolys])
        plt.legend(["currentLayerPoly","prevLayerPoly"])
        plt.axis('square')
        plt.show()
        warnings.warn(f"Layer {self.layernumber}: No intersection with prevLayer External Perimeter detected")
        return None,None

    def mergePolys(self,thesepolys:list=None)-> list:
        if not thesepolys:
            thesepolys=self.polys
        mergedPolys = unary_union(thesepolys)
        #print("Merged Geometry Type:",mergedPolys.geom_type)
        if mergedPolys.geom_type=="Polygon":
            thesepolys=[mergedPolys]
        elif mergedPolys.geom_type=="MultiPolygon" or mergedPolys.geom_type=="GeometryCollection":
            thesepolys=[poly for poly in mergedPolys.geoms]
        return thesepolys

    def spotFeaturePoints(self,featureName:str,splitAtWipe=False,includeRealStartPt=False, splitAtTravel=False)->list:
        parts=[]
        for idf,fe in enumerate(self.features):
            ftype=fe[0]
            lines=fe[1]
            start=fe[2]
            pts=[]
            isWipeMove=False
            travelstr=f"F{self.parameters.get('travel_speed')*60}"
            if featureName in ftype:
                if includeRealStartPt and idf>0:
                    sp=self.getRealFeatureStartPoint(idf)
                    if sp:pts.append(sp)
                for line in lines:
                    if "G1" in line and (not isWipeMove):
                        if (not "E" in line) and travelstr in line and splitAtTravel:
                            #print(f"Layer {self.layernumber}: try to split feature. No. of pts before:",len(pts))
                            if len(pts)>=2:#make at least 1 ls
                                parts.append(pts)
                            pts=[]# update self.features... TODO
                        elif "E" in line:     #maybe fix error of included travel moves?
                            p=getPtfromCmd(line)
                            if p:
                                pts.append(p)
                    if 'WIPE_START' in line:
                        isWipeMove=True
                        if splitAtWipe:
                            parts.append(pts)
                            pts=[]
                    if 'WIPE_END' in line:
                        isWipeMove=False
                if len(pts)>1:#fetch last one
                    parts.append(pts)
        return parts
    def spotSolidInfill(self)->None:
        parts=self.spotFeaturePoints("Solid infill",splitAtTravel=True)
        for infillpts in parts:
            if self.verifySolidInfillPts(infillpts):
                self.sinfills.append(LineString(infillpts))
    def makePolysFromSolidInfill(self,extend:float=1)->None:
        self.solidPolys=[]
        for sInfill in self.sinfills:
            infillPoly=sInfill.buffer(extend)
            self.solidPolys.append(infillPoly)
            if self.parameters.get("plotDetectedSolidInfillPoly"):
                plot_geometry(infillPoly)
                plot_geometry(sInfill,"g")
                plt.axis('square')
                plt.show()
    def verifySolidInfillPts(self,infillpts:list)->bool:
        '''Verify SollidInfillPts by checking if >=1 of the Points is inside the desired polygon-locations.'''
        for p in infillpts:
                for poly in self.oldpolys:
                    if poly.contains(p):
                        return True
        return False

    def spotBridgeInfill(self)->None:
        parts=self.spotFeaturePoints("Bridge infill",splitAtTravel=True)
        for idf,infillpts in enumerate(parts):
            self.binfills.append(BridgeInfill(infillpts))
            
    def makePolysFromBridgeInfill(self, extend: float = 1) -> None:
        """
        Create polygonal representations from bridge infill data points by extending the geometries.
    
        This function takes line segments representing infill (typically used in areas of a print that span gaps,
        like bridges) and converts them into polygons by buffering the lines. These polygons can then be used
        for further processing, such as part of the path planning for the print head.
    
        Args:
        extend (float): The distance by which the line string will be buffered to create a polygon, providing
                        a margin around the actual line segment.
    
        """
    
        # Iterate over all bridge infills in the current layer or section
        for bInfill in self.binfills:
            infillPts = bInfill.pts  # Retrieve the points defining the infill
            infillLS = LineString(infillPts)  # Create a LineString from the infill points
            infillPoly = infillLS.buffer(extend)  # Buffer the LineString to create a polygon
    
            self.polys.append(infillPoly)  # Add the new polygon to the list of polygons
            self.associatedIDs.append(bInfill.id)  # Associate the infill ID with the polygon
    
            # Plot the created polygon if the corresponding parameter is set
            if self.parameters.get("plotDetectedInfillPoly"):
                plot_geometry(infillPoly)  # Plot the buffered polygon
                plot_geometry(infillLS, "g")  # Plot the original LineString for reference
                plt.axis('square')  # Set the plot display to a square aspect ratio
                plt.show()  # Display the plot

    def getOverhangPerimeterLineStrings(self):
        parts=self.spotFeaturePoints("Overhang perimeter",includeRealStartPt=True)
        if not parts:
            parts=self.spotFeaturePoints("Overhang wall",includeRealStartPt=True)
        if parts:
            return [LineString(pts) for pts in parts]
        else:
            return []

    def verifyinfillpolys(self, minDistForValidation: float = 0.5) -> None:
       """
       Verifies each polygon in the current layer based on proximity to overhang perimeters and other criteria.
    
       This function checks if each polygon in a layer is within a specified distance to overhang perimeters, and if so,
       validates them based on additional criteria like minimum area and bridge length.
    
       Args:
       minDistForValidation (float): The maximum distance from overhang perimeters that a polygon can be to still be considered valid.
    
       Raises:
       ValueError: If no allowed space polygon is provided to the layer object.
       """
       # Retrieve line strings that represent overhang perimeters
       overhangs = self.getOverhangPerimeterLineStrings()
    
       # Debug print the count of overhangs if enabled
       if len(overhangs) > 0 and self.parameters.get("PrintDebugVerification"):
           print(f"Layer {self.layernumber}: {len(overhangs)} Overhangs found")
       
       # Check if the allowed space for arcs has been defined
       allowedSpacePolygon = self.parameters.get("AllowedSpaceForArcs")
       if not allowedSpacePolygon:
           input(f"Layer {self.layernumber}: no allowed space Polygon provided to layer obj, unable to run script. Press Enter.")
           raise ValueError(f"Layer {self.layernumber}: no allowed space Polygon provided to layer obj")
       
       # Debug print the count of polygons if enabled
       if self.parameters.get("PrintDebugVerification"):
           print("No of Polys:", len(self.polys))
       
       # Iterate through each polygon to verify validity based on specified conditions
       for idp, poly in enumerate(self.polys):
           # Check if polygon is geometrically valid
           if not poly.is_valid and self.parameters.get("PrintDebugVerification"):
               print(f"Layer {self.layernumber}: Poly{idp} is (shapely-)invalid")
               continue
    
           # Check if the polygon is within the allowed space
           if not allowedSpacePolygon.contains(poly) and self.parameters.get("CheckForAllowedSpace"):
               if self.parameters.get("PrintDebugVerification"):
                   print(f"Layer {self.layernumber}: Poly{idp} is not in allowedSpacePolygon")
               continue
    
           # Check if the polygon meets the minimum area requirement
           if poly.area < self.parameters.get("MinArea"):
               if self.parameters.get("PrintDebugVerification"):
                   print(f"Layer {self.layernumber}: Poly{idp} has too little area: {poly.area:.2f}")
                   print(f'{poly = }')
               continue
    
           # Evaluate the polygon's distance to overhang perimeters and validate based on min bridge length
           for ohp in overhangs:
               if poly.distance(ohp) < minDistForValidation:
                   if ohp.length > self.parameters.get("MinBridgeLength"):
                       self.validpolys.append(poly)  # Add polygon to valid list if it meets criteria
                       self.deleteTheseInfills.append(idp)  # Schedule the infill for deletion
                       break
                   else:
                       if self.parameters.get("PrintDebugVerification"):
                           print(f'{ohp.length = } < {self.parameters.get("MinBridgeLength") = }')
               else:
                   if self.parameters.get("PrintDebugVerification"):
                       print(f'{poly.distance(ohp) = } > {minDistForValidation}')
                       print(f"Layer {self.layernumber}: Poly{idp} is not close enough to overhang perimeters")


    def prepareDeletion(self,featurename:str="Bridge",polys:list=None)->None:
        if not polys:
            polys=self.validpolys
        for idf,fe in enumerate(self.features):
            ftype=fe[0]
            lines=fe[1]
            start=fe[2]
            deleteThis=False
            if featurename in ftype:
                for poly in polys:
                    for line in lines:
                        p=getPtfromCmd(line)
                        if p:
                            if poly.contains(p):
                                deleteThis=True
                                break
                        if deleteThis:
                            break
                if deleteThis:
                    if idf<len(self.features)-1:
                        end=self.features[idf+1][2]-1 # TODO: prevent deletion of last travel move.
                    else:
                        end=len(self.lines)
                    self.deletelines.append([start,end])
    def exportThisLine(self,linenumber:int)->bool:
        export=True
        if len(self.deletelines)>0:
            for d in self.deletelines:
                if linenumber>=d[0] and linenumber<=d[1]:
                    export=False
        return export

    def createHilbertCurveInPoly(self,poly:Polygon):
        print("making hilbert surface")
        dimensions=2
        w=self.parameters.get("solid_infill_extrusion_width")
        a=self.parameters.get("HilbertFillingPercentage")/100
        mmBetweenTravels=(self.parameters.get("aboveArcsInfillPrintSpeed")/60)*self.parameters.get("HilbertTravelEveryNSeconds")
        minX, minY, maxX, maxY=poly.bounds
        lx=maxX-minX
        ly=maxY-minY
        l=max(lx,ly)
        #Iterationcount Math explained:
        # startpoint: l/w=number of needed segments. segments=(2**iterationcount)-1. Solved for iterationcount.
        iterationCount=int(np.ceil(np.log((a*l+w)/w)/np.log(2))) # + applied ceiling function to ensucre full coverage.
        scale=w/a#l/(2**iterationCount-1)/a
        maxidx=int(2**(dimensions* iterationCount) - 1)
        locs = decode(np.arange(maxidx), 2, iterationCount)# hilbertidx->(x,y) first argument: idx, second: dimensions, third: bits per dim
        #move the curve 1 point in the smaller direction every second layer=>web the curves in z together by overlapping.
        movX=self.layernumber%2*w/a
        movY=self.layernumber%2*w/a
        x=locs[:,0]*scale+minX-movX
        y=locs[:,1]*scale+minY-movY
        hilbertPointsRaw=[[xi,yi] for xi,yi in zip(x.tolist(),y.tolist())]
        noEl=int(np.ceil(mmBetweenTravels/scale))
        buff=[]
        compositeList=[]
        #divide in subset of n elements and shuffle them to prevent localized overheating.
        for el in hilbertPointsRaw:
            p=Point(el)
            if p.within(poly):
                buff.append(p)
            else:
                if len(buff)>5:#neglegt very small pieces
                    if len(buff)>noEl*1.7:
                        compositeList.extend([buff[x:x+noEl] for x in range(0, len(buff),noEl)])
                    else:
                        compositeList.append(buff)
                buff=[]#delete single pts if there.
        if len(buff)>5:
            compositeList.append(buff) #catch last one
        random.shuffle(compositeList)
        return compositeList
    def isClose2Bridging(self,line:str,minDetectionDistance:float=3):
        if not "G1" in line:
            return False
        p=getPtfromCmd(line)
        if not p:
            return False
        if not self.lastP:
            self.lastP=Point(p.x-0.01,p.y-0.01)
        ls=LineString([p,self.lastP])
        self.lastP=p
        for poly in self.oldpolys:
            if ls.distance(poly)<minDetectionDistance:
                return True
        return False
    def spotFanSetting(self,lastfansetting:float):
        for line in self.lines:
            if "M106" in line.split(";")[0]:
                svalue=line.strip("\n").split(";")[0].split(" ")[1]
                self.fansetting=float(svalue[1:])
                return self.fansetting
        self.fansetting=lastfansetting
        return lastfansetting




class Arc():
    def __init__(self,center:Point,r:float,kwargs:dict={}) -> None:
        self.center=center
        self.r=r
        self.pointsPerCircle=kwargs.get("PointsPerCircle",80)
        self.parameters=kwargs
    def setPoly(self,poly:Polygon)->None:
        self.poly=poly
    def extractArcBoundary(self):
        circ=create_circle(self.center,self.r,self.pointsPerCircle)
        trueArc=self.poly.boundary.intersection(circ.boundary.buffer(1e-2))
        if trueArc.geom_type=='MultiLineString':
            merged=linemerge(trueArc)
        elif trueArc.geom_type=='LineString':
            self.arcline=trueArc
            return trueArc
        else:
            #print("Other Geom-Type:",trueArc.geom_type)
            merged=linemerge(MultiLineString([l for l in trueArc.geoms if l.geom_type=='LineString']))
        if merged.geom_type=="LineString":
            self.arcline=merged
            return merged
        elif merged.geom_type=="MultiLineString":
            arcList=[]
            for ls in merged.geoms:
                arc=Arc(self.center,self.r,self.parameters)
                arc.arcline=ls
                arcList.append(arc)
            return arcList
        else:
            input("ArcBoundary merging Error.Unable to run script. Press Enter.")
            raise ValueError("ArcBoundary merging Error")
    def generateConcentricArc(self,startpt:Point,remainingSpace:Polygon)->Polygon:
        circ=create_circle(startpt,self.r,self.pointsPerCircle)
        arc=circ.intersection(remainingSpace)
        self.poly=arc
        return arc

class BridgeInfill():
    def __init__(self,pts=[],id=random.randint(1,int(1e10))) -> None:
        self.pts=pts
        self.deleteLater=False
        self.id=id
