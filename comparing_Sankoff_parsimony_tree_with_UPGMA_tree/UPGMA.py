import BinaryTree as BT
class UPGMA:
    def __init__(self, dist_matrix, names = None):
        self.dist_matrix = dist_matrix
        self.leaves = []
        for i in range(len(self.dist_matrix)):
            t = BT.BinaryTree(i)
            t.distance = 0 # distance is the same as height in pseudocode
            if names: t.char = names[i]
            self.leaves.append(t)

    def get(self, i, j):
        '''get element from a lower trangular matrix'''
        return self.dist_matrix[i][j] if i > j else self.dist_matrix[j][i]

    def find_min_indexes(self):
        a, b = 1, 0
        m = self.dist_matrix[a][b]
        for i in range(1, len(self.dist_matrix)):
            for j in range(i): # weird indexing because of the lower trangular matrix
                if self.dist_matrix[i][j] < m: 
                    m = self.dist_matrix[i][j]
                    a, b = i, j
        return a, b
                
    def	cluster(self):
        '''
        Build a phylogenetic tree with evolutionary distance
        by repeatitvely coagulating the minimum distance node 
        in a distance matrix.
        Return the root node of the resulted phylogentic tree.
        '''
        while len(self.dist_matrix) > 1:
            i, j = self.find_min_indexes()
            li, lj = self.leaves.pop(i), self.leaves.pop(j)
            parent = BT.BinaryTree(-1, self.get(i,j)/2.0, li, lj)
            self.leaves.append(parent)
            
            dists = []
            for x in range(len(self.dist_matrix)):
                if x != i and x != j:
                    nc1 = li.get_num_child()
                    nc2 = lj.get_num_child()
                    dists.append((nc1*self.get(i,x)+nc2*self.get(j,x))/(nc1 + nc2))

            del self.dist_matrix[i] # remove i,j row
            del self.dist_matrix[j]
            for row in self.dist_matrix:
                del row[i] # romove i, j column
                del row[j]

            self.dist_matrix.append(dists) 
            [row.append(0.0) for row in self.dist_matrix]        
        return parent


if __name__ == '__main__':
    m = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [27, 31, 0.0, 0.0, 0.0, 0.0, 0.0],
    [8, 18, 26, 0.0, 0.0, 0.0, 0.0],
    [33, 36, 41, 31, 0.0, 0.0, 0.0],
    [18, 1, 32, 17, 35, 0.0, 0.0],
    [13, 13, 29, 14, 28, 12, 0.0]]

    matrix = UPGMA(m)
    print('*' * 60)
    print('The Phylogenetic tree is:\n')
    matrix.cluster().print_tree()

'''
Root - Dist.: 21.625
    Left - Dist.: 11.5
        Left- value:5 -- macaw
        Right- value:0 -- chiken
    Right - Dist.: 10.666666666666666
        Left - Dist.: 2.0
            Left - Dist.: 1.0
                Left- value:4 -- crab-eating monkey
                Right- value:2 -- monkey
            Right- value:1 -- human
        Right- value:3 -- cat
'''