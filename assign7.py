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
compute_midpoint()
# (3,8)
# (4,8)
# (5,7)
# (6,7)
# (7,6)
# (8,6)
# (9,5)
