import time
import tracemalloc
import sys

def main(argv):

    tracemalloc.start()
    sequences = {}

    fname = argv[1]
    prev_seq = None
    with open(fname) as f:

        start = time.time()
        count = 0
        valid_lines = 0
        lines = None
        i = 0
        while (len(line := f.readline()) > 0):
            if (line[0] == '>'):
                header = line[1:-1]
                sequences[header] = {
                    "id": header,
                    "seq": "",
                    "len": 0
                }
                # lines = [None]*1022577
                lines = []
                prev_seq = header
            elif prev_seq != None:
                if line[0] != "N":
                    valid_lines += 1
                    i += 1
                    # lines[i] = line[:-1]
                lines.append(line[:-1])
            count += 1
    sequences[header]["seq"] = "".join(lines)
    sequences[header]["len"] = len(sequences[header]["seq"])

    end = time.time()
    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')

    print("Content Lines:", valid_lines)
    print("List Size:", len(lines))
    print("Total:", valid_lines*80, "bases")

    print("Number of sequences:", len(sequences))


    print(f"Num of sequences: {len(sequences)}")
    print("Sequence Length:", sequences[header]["len"])

    print(f"Elapsed time: {(end - start) :4f}s")

    print("[ Top 10 ]")
    for stat in top_stats[:10]:
        print(stat)

if __name__ == "__main__":
    main(sys.argv)