def compute_midpoint():
    x1, y1 = 3, 8
    x2, y2 = 9, 5

    x, y = x1, y1

    dx = x2 - x1
    dy = abs(y2 - y1)

    f = -2 * dy + dx
    for _ in range(dx + 1):
        print(f"({x},{y})")
        x += 1
        if f < 0:
            y -= 1
            f += 2 * dx
        f -= 2 * dy


def rasterize_triangle(p1, p2, p3):
    # Sort vertices in counter-clockwise order
    vertices = sort_vertices(p1, p2, p3)

    # Divide mesh into two subtriangles
    subtriangles = divide_mesh(vertices)

    for triangle in subtriangles:
        # Construct bounding box for the subtriangle
        bounding_box = calculate_bounding_box(triangle)

        # Iterate pixels from left to right, bottom to top
        for x in range(bounding_box.left, bounding_box.right + 1):
            for y in range(bounding_box.bottom, bounding_box.top + 1):
                # Calculate barycentric coordinates
                barycentric_coords = calculate_barycentric_coords(triangle, x, y)

                # Check if the pixel center lies inside the triangle
                if is_inside_triangle(barycentric_coords):
                    # Color the pixel
                    color_pixel(x, y)
