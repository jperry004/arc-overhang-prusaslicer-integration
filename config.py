from shapely import Polygon

def makeFullSettingDict(gCodeSettingDict:dict) -> dict:
    """Merge Two Dictionarys and set some keys/values explicitly"""
    #the slicer-settings will be imported from GCode. But some are Arc-specific and need to be adapted by you.
    AddManualSettingsDict={
        #adapt these settings as needed for your specific geometry/printer:
        "CheckForAllowedSpace":False,# use the following x&y filter or not
        "AllowedSpaceForArcs": Polygon([[0,0],[500,0],[500,500],[0,500]]),#have control in which areas Arcs shall be generated
        "ArcCenterOffset":2, # Unit:mm, prevents very small Arcs by hiding the center in not printed section. Make 0 to get into tricky spots with smaller arcs.
        "ArcMinPrintSpeed":0.5*60,#Unit:mm/min
        "ArcPrintSpeed":1.5*60, #Unit:mm/min
        #"ArcPrintTemp":gCodeSettingDict.get("temperature"), # unit: Celsius
        "ArcTravelFeedRate":30*60, # slower travel speed, Unit:mm/min
        "ExtendIntoPerimeter":2*gCodeSettingDict.get("perimeter_extrusion_width"), #min=0.5extrusionwidth!, extends the Area for arc generation, put higher to go through small passages. Unit:mm
        #Control how much bumpiness you allow between arcs and perimeter. lower will follow perimeter better, but create a lot of very small arcs. Should be more that 1 Arcwidth! Unit:mm
        "MaxDistanceFromPerimeter":2*gCodeSettingDict.get("perimeter_extrusion_width"),
        "MinArea":35,#Unit:mm2
        "MinBridgeLength":5,#Unit:mm
        "Path2Output":"", #leave empty to overwrite the file or write to a new file. Full path required.
        "RMax":5, # the max radius of the arcs.
        "RMinStartMultiple":3, # multiple of the nozzle diameter for starting
        "TimeLapseEveryNArcs": 0, #deactivate with 0, inserts M240 after N ArcLines, 5 is a good value to start.

        #Special cooling to prevent warping:
        "aboveArcsFanSpeed":25, #0->255, 255=100%
        "aboveArcsInfillPrintSpeed":10*60, # Unit :mm/min
        "aboveArcsPerimeterFanSpeed":25, #0->255, 255=100%
        "aboveArcsPerimeterPrintSpeed":3*60, #Unit: mm/min
        "applyAboveFanSpeedToWholeLayer":True,
        "CoolingSettingDetectionDistance":5, #if the gcode line is closer than this distance to an infill polygon the cooling settings will be applied. Unit:mm
        "specialCoolingZdist":3, #use the special cooling XX mm above the arcs. Set to a negative value to disable (not recommended).

        #advanced Settings, you should not need to touch these.
        "ArcExtrusionMultiplier":0.9,
        "ArcSlowDownBelowThisDuration":3,# Arc Time below this Duration =>slow down, Unit: sec
        "ArcWidth":gCodeSettingDict.get("nozzle_diameter")*0.95, #change the spacing between the arcs,should be nozzle_diameter
        "ArcFanSpeed":255,#cooling to full blast=255
        "CornerImportanceMultiplier":0.2, # Startpoint for Arc generation is chosen close to the middle of the StartLineString and at a corner. Higher=>Cornerselection more important.
        "DistanceBetweenPointsOnStartLine":0.1,#used for redestribution, if start fails.
        "GCodeArcPtMinDist":0.1, # min Distance between points on the Arcs to for seperate GCode Command. Unit:mm
        "ExtendArcDist":1.0, # extend Arcs tangentially for better bonding bewteen them, only end-piece affected(yet), Unit:mm
        "HilbertFillingPercentage":100, # infillpercentage of the massive layers with special cooling. Uses Hilbert Curve, works not quite right yet.
        "HilbertInfillExtrusionMultiplier":1.05,
        "HilbertTravelEveryNSeconds":6, # when N seconds are driven it will continue printing somewhere else (very rough approx).
        "MinStartArcs":2, # how many arcs shall be generated in first step
        "PointsPerCircle":800, # each Arc starts as a discretized circle. Higher will slow down the code but give more accurate results for the arc-endings.
        "SafetyBreak_MaxArcNumber":2000, #max Number of Arc Start Points. prevents While loop form running for ever.
        "WarnBelowThisFillingPercentage":90, # fill the overhang at least XX%, else send a warning. Easier detection of errors in small/delicate areas. Unit:Percent
        "UseLeastAmountOfCenterPoints":True, # always generates arcs until rMax is reached, divide the arcs into pieces in needed. reduces the amount of centerpoints.

        #settings for easier debugging:
        "plotStart":False, # plot the detected geoemtry in the prev Layer and the StartLine for Arc-Generation, use for debugging
        "plotArcsEachStep":False, #plot arcs for every filled polygon. use for debugging
        "plotArcsFinal":True, #plot arcs for every filled polygon, when completely filled. use for debugging
        "plotDetectedInfillPoly":False, # plot each detected overhang polygon, use for debugging.
        "plotEachHilbert":False,
        "PrintDebugVerification":True
        }
    print(f'{gCodeSettingDict.get("nozzle_diameter") = }')
    gCodeSettingDict.update(AddManualSettingsDict)
    return gCodeSettingDict