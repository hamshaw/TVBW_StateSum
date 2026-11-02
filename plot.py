
import sympy as sp
import matplotlib.pyplot as plt

# Path to your input file
filename = "output.txt"

ranks = []
sums = []

with open(filename, "r", encoding="utf-8") as file:
    for line_num, line in enumerate(file, start=1):
        if line_num %3 != 0:
            continue
        line = line.strip()

        if not line:
            continue

        # Split into three parts:
        # entry code, rank, and the remaining algebraic expression
        parts = line.split(maxsplit=2)

        if len(parts) != 3:
            print(f"Skipping malformed line {line_num}: {line}")
            continue

        entry_code, rank_str, expression = parts

        try:
            rank = int(rank_str)

            # Parse and numerically evaluate the SymPy expression
            value = sp.sympify(expression)
            value = float(sp.N(value, 15))
            if value > 10 or value < -10:
                print(entry_code)
                continue

            ranks.append(rank)
            sums.append(value)

        except Exception as e:
            print(f"Error parsing line {line_num}: {e}")
            continue

# Create the plot
plt.figure(figsize=(10, 6))
plt.scatter(ranks, sums, alpha=0.7)

plt.xlabel("Rank")
plt.ylabel("Sum")
plt.title("Sum vs. Rank")
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("plot_no_ol.jpg")
plt.show()

