

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strings must be of the same length")
    
    distance = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            distance += 1
    return distance

def main():
    # a random 50-character adn sequence
    T = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    # a random 50-character adn sequence
    P = "CTCGCTAGC"
    # compute the hamming distance between T and

    print("T:", T)
    print("P:", P)

    distances = []
    for i in range(len(T) - len(P) + 1):
        substring = T[i:i+len(P)]
        distance = hamming_distance(substring, P)
        distances.append((i, distance))

    for i, distance in distances:
        print(f"Distance at position {i}: {distance}")

    # make a heatmap of the distances like in a continuom row from beginning to end of T.
    # closer to 0 is better, closer to 10 is worse.
    # use seaborn to make a 2d, but really is a 1d heatmap where horizontally is the position in T and color is distance at position.
    

    # we have distances of 0, 1 ... len(T) - len(P) + 1
    # indices are X
    # height is 1 to 5
    # color is distance
    dataset = []
    for i, distance in distances:
        dataset.append((i, 0, distance))
        dataset.append((i, 1, distance))
        dataset.append((i, 2, distance))
        dataset.append((i, 3, distance))
        dataset.append((i, 4, distance))

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    df = pd.DataFrame(dataset, columns=["Position", "Height", "Distance"])
    plt.figure(figsize=(20, 5))
    sns.heatmap(df.pivot("Height", "Position", "Distance"), cmap="coolwarm", cbar_kws={'label': 'Hamming Distance'})
    plt.title("Hamming Distance Heatmap")
    plt.xlabel("Position in T")
    plt.ylabel("Height")
    plt.show()
    




if __name__ == "__main__":
    main()