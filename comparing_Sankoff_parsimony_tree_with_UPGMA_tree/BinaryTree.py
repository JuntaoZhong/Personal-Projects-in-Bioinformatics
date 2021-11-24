class BinaryTree:
    def __init__(self, val, distance = 0, left = None, right = None, parent = None, char = None,source = None,temp = None):
        self.value = val
        self.distance = distance
        self.left = left
        self.right = right
        self.parent = parent
        self.char = char
        self.source = source
        self.temp = temp

    def print_tree(self):
        self.__print_tree_rec(0, 'Root')

    def __print_tree_rec(self, level, side):
        to_print = '\t'*level + side
        if self.value >= 0:
            to_print += '- value:'+ str(self.value)
            if self.char: to_print += (" -- " + self.char)
            if self.distance: to_print += (" - dist: " + self.distance)
            print(to_print) 
        else:
            print(to_print, '- Dist.:', self.distance)
            if(self.left != None):
                self.left.__print_tree_rec(level+1, 'Left') 
            if(self.right != None):
                self.right.__print_tree_rec(level+1, 'Right')
    
    def print_tree_char(self):
        out = []
        self.__print_tree_char(0, 'Root', out)
        return out

    def __print_tree_char(self, level, side, out):
        min_key = min(self.distance.keys(), key=(lambda k: self.distance[k]))
        out.append(min_key)
        if self.value <= 0:
            if(self.left != None):
                self.left.__print_tree_char(level+1, 'Left', out) 
            if(self.right != None):
                self.right.__print_tree_char(level+1, 'Right', out)

    def get_num_child(self):
        if self.value >=0: 
            return 1 # val(internal nodes) < 0
        res = 0
        if self.left: res += self.left.get_num_child()
        if self.right: res += self.right.get_num_child()
        return res


if __name__ == "__main__":
    a = BinaryTree(1)
    b = BinaryTree(2)
    c = BinaryTree(3)
    d = BinaryTree(4)
    e = BinaryTree(-1, 3, b, c)
    f = BinaryTree(-1, 3, d, a)
    g = BinaryTree(-1, 6, e, f)

    g.print_tree()
    e.print_tree()

    print(g.get_leaf())
    print(e.get_leaf())
