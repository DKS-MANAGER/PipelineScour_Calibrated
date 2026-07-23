import struct

with open('constant/triSurface/cylinder.stl', 'rb') as f:
    header = f.read(80)
    num_triangles = struct.unpack('<I', f.read(4))[0]
    x_vals = []
    y_vals = []
    z_vals = []
    for _ in range(num_triangles):
        # Read normal (3 floats)
        f.read(12)
        # Read 3 vertices (each has 3 floats: x, y, z)
        for _ in range(3):
            vx, vy, vz = struct.unpack('<fff', f.read(12))
            x_vals.append(vx)
            y_vals.append(vy)
            z_vals.append(vz)
        # Read attribute byte count (2 bytes)
        f.read(2)

print("X bounds:", min(x_vals), max(x_vals))
print("Y bounds:", min(y_vals), max(y_vals))
print("Z bounds:", min(z_vals), max(z_vals))

