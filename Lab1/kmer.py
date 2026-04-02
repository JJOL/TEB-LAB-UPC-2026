import pandas as pd
import matplotlib.pyplot as plt
import sys

def main(csv_file, top=50):
    # Load CSV
    df = pd.read_csv(csv_file)

    # Sort by count descending
    # df = df.sort_values("count", ascending=False)

    print("Total k-mers:", len(df))

    # Keep only top N (otherwise unreadable for large K)
    df_top = df.head(top)

    # Plot bar chart
    plt.figure(figsize=(14,6))
    plt.bar(df_top["kmer"], df_top["count"])

    plt.title(f"Top {top} most frequent k-mers")
    plt.xlabel("k-mer")
    plt.ylabel("Frequency")

    plt.xticks(rotation=90)
    plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_kmers.py hist.csv [topN]")
        sys.exit(1)

    csv_file = sys.argv[1]
    topN = int(sys.argv[2]) if len(sys.argv) > 2 else 50

    main(csv_file, topN)
