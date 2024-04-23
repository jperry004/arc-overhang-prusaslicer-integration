from shapely.geometry import LineString, Polygon, Point
import numpy as np
import matplotlib.pyplot as plt
import re



features = [[';TYPE:Perimeter\n', [';LAYER_CHANGE\n', ';Z:4.2\n', ';HEIGHT:0.15\n', ';BEFORE_LAYER_CHANGE\n', 'G92 E0\n', ';4.2\n', '\n', '\n', 'G1 E-4.5 F2400\n', ';WIPE_START\n', 'G1 F9600\n', 'G1 X122.118 Y116.642 E-1.5\n', ';WIPE_END\n', 'G1 Z4.2 F12000\n', 'M73 P46 R8\n', ';AFTER_LAYER_CHANGE\n', ';4.2\n', '; printing object Overhang_Test.stl id:0 copy 0\n', 'G1 Z4.2\n', 'G1 X110.953 Y121.444\n', 'G1 Z4.2\n', 'G1 E6.4 F1800\n', ';WIDTH:0.719999\n', 'G1 F2400\n', 'G1 X110.944 Y121.435 E.00055\n', 'G1 X110.932 Y119.917 E.06511\n', 'G1 X110.893 Y119.538 E.01634\n', 'G1 X110.738 Y119.16 E.01752\n', 'G1 X110.624 Y119.077 E.00605\n', 'G1 X110.739 Y118.994 E.00608\n', 'G1 X110.851 Y118.779 E.0104\n', 'G1 X110.918 Y118.412 E.016\n', 'G1 X110.917 Y116.101 E.09913\n', 'G1 X110.873 Y115.803 E.01292\n', 'G1 X110.693 Y115.47 E.01624\n', 'G1 X110.873 Y115.136 E.01627\n', 'G1 X110.917 Y114.84 E.01284\n', 'G1 X110.918 Y113.453 E.05949\n', 'G1 X110.889 Y113.21 E.0105\n', 'G1 X110.664 Y112.757 E.0217\n', 'G1 X110.869 Y112.393 E.01792\n', 'G1 X110.918 Y112.088 E.01325\n', 'G1 X110.941 Y109.548 E.10895\n', 'G1 X123.463 Y109.548 E.53711\n', 'G1 X123.463 Y115.5 E.2553\n', 'G1 X123.463 Y121.452 E.2553\n', 'G1 X111.04 Y121.452 E.53287\n', 'G1 E-4.5\n', ';WIPE_START\n', 'G1 F9600\n', 'G1 X110.953 Y121.444 E-.02075\n', 'G1 X110.944 Y121.435 E-.00302\n', 'G1 X110.932 Y119.917 E-.36054\n', 'G1 X110.893 Y119.538 E-.09049\n', 'G1 X110.738 Y119.16 E-.09703\n', 'G1 X110.624 Y119.077 E-.03349\n', 'G1 X110.739 Y118.994 E-.03368\n', 'G1 X110.851 Y118.779 E-.05758\n', 'G1 X110.918 Y118.412 E-.0886\n', 'G1 X110.917 Y116.101 E-.54886\n', 'G1 X110.873 Y115.803 E-.07154\n', 'G1 X110.693 Y115.47 E-.0899\n', 'G1 X110.702 Y115.453 E-.00452\n', ';WIPE_END\n', 'G1 X110.291 Y122.14 F12000\n', 'G1 E6.4 F1800\n'], 0], [';TYPE:External perimeter\n', ['G1 F1500\n', 'G1 X110.257 Y121.455 E.02942\n', 'G1 X110.244 Y119.922 E.06576\n', 'G1 X110.222 Y119.697 E.0097\n', 'G1 X110.166 Y119.543 E.00703\n', 'G1 X109.976 Y119.415 E.00983\n', ';WIDTH:0.698431\n', 'G1 X109.547 Y119.493 E.01812\n', ';WIDTH:0.671699\n', 'G1 X109.026 Y119.624 E.02143\n', 'G1 X109.026 Y118.936 E.02744\n', 'G1 X109.461 Y118.864 E.01758\n', ';WIDTH:0.698431\n', 'G1 X109.979 Y118.732 E.02221\n', ';WIDTH:0.719999\n', 'G1 X110.171 Y118.606 E.00985\n', 'G1 X110.231 Y118.412 E.00871\n', 'G1 X110.23 Y116.101 E.09913\n', 'G1 X110.144 Y115.874 E.01041\n', 'G1 X109.96 Y115.827 E.00815\n', ';WIDTH:0.692209\n', 'G1 X109.775 Y115.78 E.00786\n', ';WIDTH:0.673897\n', 'G1 X109.691 Y115.82 E.00372\n', 'G1 X109.477 Y115.692 E.00998\n', 'G1 X109.252 Y115.438 E.01358\n', ';WIDTH:0.681838\n', 'G1 X109.057 Y114.963 E.0208\n', ';WIDTH:0.710375\n', 'G1 X109.015 Y114.54 E.01798\n', 'G1 X109.208 Y114.809 E.014\n', ';WIDTH:0.685811\n', 'G1 X109.512 Y115.046 E.01571\n', ';WIDTH:0.673897\n', 'G1 X109.837 Y115.152 E.01368\n', ';WIDTH:0.692209\n', 'G1 X109.99 Y115.109 E.00654\n', ';WIDTH:0.719999\n', 'G1 X110.144 Y115.066 E.00686\n', 'G1 X110.199 Y114.98 E.00438\n', 'G1 X110.229 Y114.839 E.00618\n', 'G1 X110.225 Y113.397 E.06185\n', 'G1 X110.145 Y113.227 E.00806\n', 'G1 X109.939 Y113.113 E.0101\n', 'G1 X109.834 Y113.113 E.0045\n', 'G1 X109.537 Y113.188 E.01314\n', ';WIDTH:0.721891\n', 'G1 X109.2 Y113.442 E.01815\n', ';WIDTH:0.726032\n', 'G1 X109.008 Y113.793 E.01731\n', 'G1 X109.018 Y113.652 E.00612\n', ';WIDTH:0.724596\n', 'G1 X109.121 Y113.2 E.02002\n', ';WIDTH:0.719999\n', 'G1 X109.334 Y112.793 E.0197\n', 'G1 X109.654 Y112.533 E.01769\n', 'G1 X110.145 Y112.305 E.02322\n', 'G1 X110.198 Y112.225 E.00412\n', 'G1 X110.23 Y112.082 E.00629\n', 'G1 X110.258 Y109.054 E.12989\n', 'G1 X110.275 Y108.86 E.00835\n', 'G1 X123.9 Y108.86 E.58443\n', 'M204 P250\n'], 65], [';TYPE:Overhang perimeter\n', [';WIDTH:0.584808\n', ';HEIGHT:0.584808\n', 'G1 F300\n', 'G1 X124.151 Y108.86 E.02803\n', 'G1 X124.151 Y115.5 E.74151\n', 'G1 X124.151 Y122.14 E.74151\n', 'G1 X123.9 Y122.14 E.02803\n', 'M204 P500\n'], 130], [';TYPE:External perimeter\n', [';WIDTH:0.719999\n', ';HEIGHT:0.15\n', 'G1 F1500\n', 'G1 X110.381 Y122.14 E.57988\n', 'G1 E-4.5 F2400\n', ';WIPE_START\n', 'G1 F9600\n', 'G1 X110.291 Y122.14 E-.02138\n', 'G1 X110.257 Y121.455 E-.16289\n', 'G1 X110.244 Y119.922 E-.3641\n', 'G1 X110.222 Y119.697 E-.05369\n', 'G1 X110.166 Y119.543 E-.03892\n', 'G1 X109.976 Y119.415 E-.05441\n', 'G1 X109.547 Y119.493 E-.10356\n', 'G1 X109.026 Y119.624 E-.12759\n', 'G1 X109.026 Y118.936 E-.1634\n', 'G1 X109.461 Y118.864 E-.10472\n', 'G1 X109.979 Y118.732 E-.12696\n', 'G1 X110.171 Y118.606 E-.05454\n', 'G1 X110.231 Y118.412 E-.04823\n', 'G1 X110.231 Y118.094 E-.07561\n', ';WIPE_END\n', 'G1 X110.511 Y118.218 F12000\n', 'G1 X110.85 Y117.451\n', 'G1 X110.901 Y116.231\n', 'G1 X110.261 Y115.741\n', 'G1 X108.795 Y114.696\n', 'G1 X108.768 Y114.499\n', 'G1 X109.011 Y114.465\n', 'G1 E6.4 F1800\n', ';WIDTH:0.715634\n', 'G1 F1500\n', 'G1 X108.974 Y114.07 E.01691\n', ';WIDTH:0.726032\n', 'G1 X109.008 Y113.793 E.01208\n', 'G1 E-4.5 F2400\n', ';WIPE_START\n', 'G1 F9600\n', 'G1 X108.974 Y114.07 E-.06628\n', 'G1 X109.011 Y114.465 E-.09422\n', ';WIPE_END\n', 'G1 E-1.3395 F2400\n', 'G1 X109.262 Y114.402 F12000\n', 'M73 P47 R8\n', 'G1 X109.331 Y114.621\n', 'G1 X109.351 Y114.671\n', 'G1 X109.489 Y114.893\n', 'G1 X109.529 Y114.934\n', 'G1 X110.214 Y115.306\n', 'G1 X110.289 Y114.934\n', 'G1 X111.429 Y115.163\n', 'G1 E6.4 F1800\n'], 139], [';TYPE:Solid infill\n', [';WIDTH:0.808901\n', 'G1 F2400\n', 'G1 X111.459 Y115.085 E.00405\n', ';WIDTH:0.842601\n', 'G1 X111.489 Y115.008 E.00418\n', ';WIDTH:0.876302\n', 'G1 X111.519 Y114.931 E.00435\n', 'G1 X111.564 Y115.017 E.00511\n', ';WIDTH:0.838394\n', 'G1 X111.659 Y115.136 E.00766\n', ';WIDTH:0.791935\n', 'G1 X111.754 Y115.254 E.00718\n', ';WIDTH:0.745476\n', 'G1 X111.849 Y115.373 E.00677\n', ';WIDTH:0.699017\n', 'G1 X111.943 Y115.491 E.00627\n', ';WIDTH:0.664543\n', 'G1 X112.047 Y115.996 E.02033\n', 'G1 X112.05 Y117.765 E.06976\n', ';WIDTH:0.663872\n', 'G1 X111.961 Y118.261 E.01985\n', ';WIDTH:0.623366\n', 'G1 X111.951 Y118.559 E.01099\n', ';WIDTH:0.666937\n', 'G1 X111.99 Y118.876 E.01264\n', ';WIDTH:0.710507\n', 'G1 X112.03 Y119.193 E.01352\n', 'G1 X112.056 Y119.782 E.02494\n', ';WIDTH:0.663019\n', 'G1 X112.064 Y120.953 E.04607\n', ';WIDTH:0.659618\n', 'G1 X111.439 Y120.953 E.02446\n', 'G1 X111.435 Y120.643 E.01213\n', ';WIDTH:0.663019\n', 'G1 X111.426 Y119.788 E.03364\n', ';WIDTH:0.710507\n', 'G1 X111.357 Y119.256 E.02269\n', ';WIDTH:0.727101\n', 'G1 X111.288 Y119.07 E.0086\n', 'G1 X111.371 Y118.872 E.0093\n', ';WIDTH:0.7137\n', 'G1 X111.376 Y118.703 E.00719\n', ';WIDTH:0.674792\n', 'G1 X111.381 Y118.535 E.00674\n', ';WIDTH:0.635883\n', 'G1 X111.386 Y118.366 E.00637\n', ';WIDTH:0.630423\n', 'G1 X111.402 Y118.054 E.01166\n', ';WIDTH:0.663872\n', 'G1 X111.418 Y117.742 E.01231\n', ';WIDTH:0.664543\n', 'G1 X111.414 Y115.997 E.06881\n', ';WIDTH:0.666014\n', 'G1 X111.413 Y115.967 E.00119\n', ';WIDTH:0.67407\n', 'G1 X111.275 Y115.47 E.02065\n', ';WIDTH:0.714517\n', 'G1 X111.324 Y115.362 E.00505\n', ';WIDTH:0.754963\n', 'G1 X111.373 Y115.254 E.00535\n', ';WIDTH:0.795409\n', 'G1 X111.422 Y115.146 E.00564\n', ';WIDTH:0.835856\n', 'G1 X111.47 Y115.038 E.00592\n', ';WIDTH:0.876302\n', 'G1 X111.519 Y114.931 E.0062\n', ';WIDTH:0.911864\n', 'G1 X111.538 Y114.682 E.0137\n', ';WIDTH:0.947426\n', 'G1 X111.557 Y114.434 E.0142\n', ';WIDTH:0.982988\n', 'G1 X111.576 Y114.186 E.01475\n', 'G1 X111.558 Y113.898 E.01711\n', ';WIDTH:0.94475\n', 'G1 X111.54 Y113.61 E.01642\n', ';WIDTH:0.912192\n', 'G1 X111.543 Y113.461 E.00818\n', ';WIDTH:0.719999\n', 'G1 X111.254 Y112.766 E.03229\n', ';WIDTH:0.744986\n', 'G1 X111.347 Y112.957 E.00944\n', ';WIDTH:0.791018\n', 'G1 X111.439 Y113.147 E.00999\n', ';WIDTH:0.83141\n', 'G1 X111.474 Y113.252 E.00552\n', ';WIDTH:0.871801\n', 'G1 X111.509 Y113.357 E.0058\n', ';WIDTH:0.912192\n', 'G1 X111.543 Y113.461 E.006\n', 'G1 X111.62 Y113.354 E.00723\n', ';WIDTH:0.869545\n', 'G1 X111.696 Y113.247 E.00685\n', ';WIDTH:0.826897\n', 'G1 X111.773 Y113.14 E.00653\n', ';WIDTH:0.784249\n', 'G1 X111.849 Y113.032 E.00619\n', ';WIDTH:0.741602\n', 'G1 X111.925 Y112.925 E.00581\n', ';WIDTH:0.698954\n', 'G1 X112.02 Y112.607 E.0138\n', ';WIDTH:0.690372\n', 'G1 X112.052 Y112.091 E.02122\n', ';WIDTH:0.664817\n', 'G1 X112.06 Y110.362 E.06821\n', ';WIDTH:0.659437\n', 'G1 X112.061 Y110.047 E.01232\n', 'G1 X111.436 Y110.047 E.02445\n', 'G1 X111.433 Y110.358 E.01217\n', ';WIDTH:0.664817\n', 'G1 X111.42 Y112.087 E.06822\n', ';WIDTH:0.690643\n', 'G1 X111.37 Y112.448 E.01497\n', 'G1 E-6\n', 'G1 X111.37 Y112.448 F12000\n', 'G1 X112.707 Y118.978\n', 'G1 E6.4 F1800\n'], 183], [';TYPE:Internal infill\n', [';WIDTH:0.72\n', 'G1 F3000\n', 'G1 X112.581 Y118.418 E.02462\n', 'G1 X112.709 Y117.809 E.02669\n', 'G1 X112.71 Y115.87 E.08317\n', 'G1 X112.643 Y115.484 E.0168\n', 'G1 X112.515 Y115.166 E.0147\n', 'G1 X113.878 Y110.077 E.22598\n', 'G1 X114.277 Y110.077 E.01711\n', 'G1 X121.878 Y117.678 E.46108\n', 'G1 X121.878 Y116.866 E.03483\n', 'G1 X112.713 Y119.321 E.40698\n', 'G1 X112.722 Y120.923 E.06872\n', 'G1 X115.247 Y120.923 E.10831\n', 'G1 E-6 F2400\n', 'G1 X115.247 Y120.923 F12000\n', 'G1 X122.486 Y121.003\n', 'G1 E6.4 F1800\n'], 307], [';TYPE:Solid infill\n', [';WIDTH:0.559912\n', 'G1 F2400\n', 'G1 X123.013 Y121.003 E.01734\n', 'G1 X123.013 Y109.997 E.36221\n', 'G1 X122.486 Y109.997 E.01734\n', 'G1 X122.486 Y120.659 E.35089\n', 'M73 P47 R7\n', '; stop printing object Overhang_Test.stl id:0 copy 0\n'], 326]]


def create_buffer(line, distance, side='left'):
    # Assuming 'left' or 'right' side buffer based on a simple 2D Cartesian coordinate system
    original_coords = list(line.coords)
    offset_coords = []

    for i in range(len(original_coords) - 1):
        p1 = np.array(original_coords[i])
        p2 = np.array(original_coords[i + 1])
        
        # Calculate the direction vector from p1 to p2
        direction = p2 - p1
        # Calculate the normal vector (perpendicular)
        if side == 'left':
            normal = np.array([-direction[1], direction[0]])
        elif side == 'right':
            normal = np.array([direction[1], -direction[0]])
        
        # Normalize the normal vector
        normal = normal / np.linalg.norm(normal)
        
        # Create new points offset by the distance in the direction of the normal
        new_p1 = p1 + normal * distance
        new_p2 = p2 + normal * distance
        
        offset_coords.append(tuple(new_p1))
        if i == len(original_coords) - 2:  # Also append the last point
            offset_coords.append(tuple(new_p2))

    # Construct the polygon that includes original and offset coordinates
    # The path will trace along the original coordinates and back along the offset ones
    if side == 'left':
        all_coords = original_coords + offset_coords[::-1]  # Reverse offset coordinates to loop back correctly
    else:
        all_coords = offset_coords + original_coords[::-1]

    if all_coords[0] != all_coords[-1]:
        all_coords.append(all_coords[0])  # Ensure the polygon is closed

    return Polygon(all_coords)


def plot_polygon(polygon):
    # Extract exterior coordinates of the polygon
    x, y = polygon.exterior.xy
    
    fig, ax = plt.subplots()
    ax.fill(x, y, alpha=0.5, color='orange', edgecolor='black')  # Fill the polygon with a color and edge
    
    # Set plot limits to the bounds of your polygon, with some padding
    x_min, y_min, x_max, y_max = polygon.bounds
    ax.set_xlim(x_min - 1, x_max + 1)
    ax.set_ylim(y_min - 1, y_max + 1)
    
    # Set aspect of the plot to be equal, so the polygon isn't skewed
    ax.set_aspect('equal')
    
    plt.grid(True)
    plt.title('Shapely Polygon Plot')
    plt.show()

# # Example usage
# line = LineString([(0, 0), (2, 2)])
# buffered_polygon = create_buffer(line, 1, side='right')
# print(buffered_polygon)
# plot_polygon(buffered_polygon)


def parse_gcode(gcode_lines):
    coords = []
    # Updated regex pattern to capture X and Y values flexibly.
    # This pattern is adjusted to handle optional spaces and any characters that may appear between 'G1' and the coordinates.
    pattern = re.compile(r'G1\s+[^XY]*X([-\d.]+)\s+[^XY]*Y([-\d.]+)', re.IGNORECASE)
    
    for line in gcode_lines:
        # print("Processing line:", line.strip())  # Ensure we see what's being processed
        
        match = pattern.search(line)
        if match:
            x = float(match.group(1))
            y = float(match.group(2))
            coords.append((x, y))
            # print(f"Match found: X={x}, Y={y}")  # Debugging to confirm matches    
    return coords

def decide_buffer_direction(feature_coords, all_coords, distance):
    line = LineString(feature_coords)
    left_buffer = line.buffer(distance, single_sided=True)
    right_buffer = line.buffer(-distance, single_sided=True)

    left_intersect = sum(1 for point in all_coords if left_buffer.contains(Point(point)))
    right_intersect = sum(1 for point in all_coords if right_buffer.contains(Point(point)))

    if left_intersect > right_intersect:
        return 'left'
    elif right_intersect > left_intersect:
        return 'right'
    else:
        return 'undecided'
all_feature_coords = []
for feature in features:
    # Ensure feature[1] (assumed to be a list of G-code lines) is passed as a whole to parse_gcode
    feature_coords = parse_gcode(feature[1])
    # print(len(feature_coords))
    # if len(feature_coords) < 10:
    #     print(feature_coords)
    all_feature_coords.extend(feature_coords)

# print(all_feature_coords)


# Analyze each feature and decide buffer direction
buffer_decisions = []
for feature in features:
    feature_coords = parse_gcode(feature[1])
    # print(feature_coords)
    buffer_side = decide_buffer_direction(feature_coords, all_feature_coords, 0.5)
    buffer_decisions.append((feature[0], buffer_side))

# Output decisions
for decision in buffer_decisions:
    print("Feature Type:", decision[0], "Buffer Direction:", decision[1])








