from graphviz import Digraph

calls_count = 0
call_id = 0
graph = None
current_parent = None

def main():
    global graph, calls_count, call_id, current_parent
    
    print("Hello, World!")
    x = "ATGATGAAGTGTGTAGCTCGCGGCCGATCGACTGCACGTACGTAGCCGCGACGATCTAGCTATATATAGCTAGTCGATCGCGACTGCATGCATCGCCGTCGCTCCTCCCGAATAACTAGCTACAGATAGAGAGAGAGATCGACTAGCTACGATCGCACTGTTTGACCACTGCAG"
    y = "AGCTAGTCGATCGCGACTGCATGCATCGCCGTCGCTCCTCCCGAATAACTAGCTACAGATAGAGAGAGAGATCGACTAGCTACGATCGCACTGTTTGACCACTGCAG"
    # x = "Juan Jose Olivera"
    # y = "JuanJose"
    # x = "Juan Jose"
    # y = "JuanJose"
    x = "ACT"
    y = "AGT"
    
    # Initialize graph
    # graph = Digraph(comment='Edit Distance Recursive Calls')
    # graph.attr(rankdir='TB')
    # graph.attr('node', shape='box', style='rounded,filled', fillcolor='lightblue')
    
    calls_count = 0
    call_id = 0
    current_parent = None
    
    result = dp_edit_distance(x, y)
    
    print(f"Edit distance between '{x}' and '{y}': {result}")
    print(f"Number of calls to edit_distance: {calls_count}")
    
    # Save and display graph
    # output_file = 'edit_distance_calls'
    # graph.render(output_file, format='pdf', cleanup=True)
    # print(f"Graph saved to {output_file}.pdf")

def edit_distance(s1: str, s2: str, parent_id=None) -> int:
    global calls_count, call_id, graph, current_parent
    
    calls_count += 1
    my_id = call_id
    call_id += 1
    
    # Create node for this call
    label = f"ed('{s1}', '{s2}')"
    graph.node(str(my_id), label)
    
    # Add edge from parent if exists
    if parent_id is not None:
        graph.edge(str(parent_id), str(my_id))
    
    m = len(s1)
    n = len(s2)

    if m == 0:
        return n
    if n == 0:
        return m
    
    cost = edit_distance(s1[1:], s2[1:], my_id)
    if s1[0] != s2[0]:
        cost += 1

    insert_cost = edit_distance(s1, s2[1:], my_id) + 1
    delete_cost = edit_distance(s1[1:], s2, my_id) + 1
    return min(cost, insert_cost, delete_cost)


def dp_edit_distance(s1: str, s2: str) -> int:
    m = len(s1)
    n = len(s2)
    
    dp = [[(0,-1)] * (n + 1) for _ in range(m + 1)]
    
    for i in range(m + 1):
        dp[i][0] = (i, -1)
    for j in range(n + 1):
        dp[0][j] = (j, -1)
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # indicate the operation (cell origin) that led to the minimum cost to trace back the path
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = (dp[i - 1][j - 1][0], 0)  # no cost, from diagonal
            else:
                dp[i][j] = (dp[i - 1][j - 1][0] + 1, 0)  # substitution, from diagonal

            if dp[i][j][0] > dp[i][j - 1][0] + 1:
                dp[i][j] = (dp[i][j - 1][0] + 1, 1)  # insertion, from left
            if dp[i][j][0] > dp[i - 1][j][0] + 1:
                dp[i][j] = (dp[i - 1][j][0] + 1, 2)  # deletion, from top

    # reconstruct path (optional, for visualization)
    path = []
    i, j = m, n
    while i > 0 or j > 0:
        path.append((i, j))
        if dp[i][j][1] == 0:
            i -= 1
            j -= 1
        elif dp[i][j][1] == 1:
            j -= 1
        else:
            i -= 1

    path.reverse()  # optional, to have the path from start to end

    return dp[m][n][0], path

if __name__ == "__main__":
    main()