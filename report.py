import matplotlib.pyplot as plt


# Times collected from running the program
data = {
    "armadillo_small": {"size": 3112, "time": 12.2288},
    "lucy_small": {"size": 2804, "time": 13.9946},
    "bunny_small": {"size": 1392, "time": 5.54895},
    "armadillo": {"size": 345944, "time": 17.043},
    "lucy": {"size": 2805572, "time": 33.7117},
    "bunny": {"size": 69451, "time": 8.2359},
}

names = list(data.keys())
sizes = [entry["size"] for entry in data.values()]
times = [entry["time"] for entry in data.values()]

plt.figure(figsize=(10, 6))
plt.scatter(sizes, times, marker="o", color="b")

for i, name in enumerate(names):
    plt.text(sizes[i] + 100, times[i], name, fontsize=9)

plt.title("Size of Mesh vs Time to Render")
plt.xlabel("Size")
plt.ylabel("Time")
plt.grid(True)
plt.show()
