#To install tqdm run sage -pip install tqdm in a terminal
from tqdm import tqdm

#===============================================================================
#        Fonction à ajouter à Implicit Suffix Tree
#===============================================================================
def LZ_decomposition(self):
    r"""
    Return the Lempel-Ziv decomposition of the self.word() in the form of a list
    iB of index such that the blocks of the decomposition are 
    self.word()[iB[k]:iB[k+1]]

    The Lempel-Ziv decomposition is the factorisation u_1...u_k of a word 
    w=x_1...x_n such that u_i is the longest prefix of u_i...u_k that has an
    occurence starting before u_i or a letter if the prefix is empty.
    
    EXAMPLE:

        sage: w = Word('abababb')
        sage: T = w.suffix_tree()
        sage: T.LZ_decomposition()
        [0, 1, 2, 6, 7]
        sage: w = Word('abaababacabba')
        sage: T = w.suffix_tree()
        sage: T.LZ_decomposition()
        [0, 1, 2, 3, 6, 8, 9, 11, 13]
        sage: w = Word([0,0,0,1,1,0,1])
        sage: T = w.suffix_tree()
        sage: T.LZ_decomposition()
        [0, 1, 3, 4, 5, 7]
        sage: w=Word('0000100101')
        sage: T=w.suffix_tree()
        sage: T.LZ_decomposition()
        [0, 1, 4, 5, 9, 10]
    """
    iB=[0]
    i=0
    w=self.word()
    while i<len(w):
        l=0
        ((x,y),successor)=self._find_transition(0, w[i])
        x=x-1
        while x<i+l:
            if y==None:
                l=len(w)-i
            else:
                l+=y-x
            if i+l>=len(w):
                l=len(w)-i
                break
            ((x,y),successor)=self._find_transition(successor,w[i+l])
            x=x-1
        i+=max(1,l)
        iB.append(i)
    return iB

def leftmost_covering_set(self):
    r"""
    Compute the leftmost covering set of squares pair in self.word(). Return 
    square as pair (i,l) specifying self.word()[i:i+l]

    A leftmost covering set is a set such that the leftmost occurence (j,l) of a
    type of square in self.word() is covered by a pair (i,l) in the set for all
    types of squares. We say that (j,l) is covered by (i,l) if (i,l), (i+1,l),
    ..., (j,l) are all squares.

    The set is return in the form of a list P such that P[i] contains all the
    the length of square starting at i in the set. The list P[i] are sort in
    decreasing order.
    
    EXAMPLES:

        sage: w=Word('abaabaabbaaabaaba')
        sage: T=w.suffix_tree()
        sage: T.leftmost_covering_set()
        [[6], [6], [2], [], [], [], [], [2], [], [], [6, 2], [], [], [], [], [], []]
        sage: w=Word('abaca')
        sage: T=w.suffix_tree()
        sage: T.leftmost_covering_set()
        [[], [], [], [], []]

    REFERENCE:

    * [1] Gusfield, D., & Stoye, J. (2004). Linear time algorithms for finding
      and representing all the tandem repeats in a string. Journal of Computer
      and System Sciences, 69(4), 525-546.
    """
    def condition1_square_pairs(i):
        r"""
        Compute the square that has their center in the i-th block of 
        LZ-decomposition and that start in the i-th block and end in the
        (i+1)-th
        """
        for k in range(1,B[i+1]-B[i]+1):
            q=B[i+1]-k
            k1=w.longest_forward_extension(B[i+1],q)
            k2=w.longest_backward_extension(B[i+1]-1,q-1)
            start=max(q-k2,q-k+1)
            if k1+k2>=k and k1>0 and start>=B[i]:
                #print "Condition 1 yield (%s,%s) for block %s" %(start,2*k,i)
                #print start>=B[i]
                yield (start,2*k)
    
    def condition2_square_pairs(i):
        r"""
        Compute the squares that has their center in the i-th block of the
        LZ-decomposition and that starts in the (i-1)-th block or before. Their
        end is either in the i-th or the (i+1)-th block
        """
        try:
            end=B[i+2]-B[i]+1
        except IndexError:
            end=B[i+1]-B[i]+1
        for k in range(2,end):
            q=B[i]+k
            k1=w.longest_forward_extension(B[i],q)
            k2=w.longest_backward_extension(B[i]-1,q-1)
            start=max(B[i]-k2,B[i]-k+1)
            #print "k=%s q=%s k1=%s k2=%s, start=%s, h1=%s" %(k,q,k1,k2,start,B[i+1])
            if k1+k2>=k and k1>0 and start+k<=B[i+1] and k2>0:# and start<B[i]:
                #print "Condition 2 yield (%s,%s) for block %s" %(start,2*k,i)
                #print start<B[i]
                yield (start,2*k)

    w=self.word()
    B=self.LZ_decomposition()
    P=[[] for _ in w]
    for i in range(len(B)-1):
        squares=list(condition2_square_pairs(i))+list(condition1_square_pairs(i))
        for (i,l) in squares:
            P[i].append(l)
    for l in P:
        l.reverse()
    return P
    
def count_and_skip(self,node,(i,j)):
    r"""
    Use count and skip trick to follow the path starting at "node" and
    reading self.word()[i:j]. We assume that reading self.word()[i:j] is 
    possible from "node"

    INPUTS:

        node - explicit node of T
        (i,j) - Indices of factor T.word()[i:j]

    OUTPUT:

        The node obtained by starting at "node" and following the edges
        labeled by the letter of T.word()[i:j]. Returns "("explicit",
        end_node) if w ends at a "end_node", and "("implicit", edge, d)" if
        it ends at a spot along an edge.

    EXAMPLES:

        sage:T=Word('00110111011').suffix_tree()
        sage:T.count_and_skip(5,(2,5))
        ("implicit", (9, 10), 2)
        sage:T.count_and_skip(0, (1, 4))
        ("explicit", 7)
    """
    if i==j: #We're done reading the factor
        return ('explicit',node)
    transition=self._find_transition(node,self._letters[i])
    child=transition[1]
    if transition[0][1]==None: #The child is a leaf
        edge_length=len(self.word())-transition[0][0]+1
    else:
        edge_length=transition[0][1]-transition[0][0]+1
    if edge_length>j-i: #The reading stop on this edge
        return ('implicit',(node,child),j-i)
    return self.count_and_skip(child,(i+edge_length,j))

def suffix_walk(self,(edge,l)):
    r"""
    Compute the suffix walk from the input state. If the input state is path
    label "aw" with "a" a letter, the output is the state of "w".

    INPUTS:

        edge - the edge containign the state
        l - the string-depth of the state on edge (l>0)

    OUTPUT:

        Returns "("explicit", end_node) if the state of "w" is an explicit state
        and "("implicit", edge, d)" if the state of "w" is implicit on "edge".

    EXAMPLES:

        sage:T=Word('00110111011').suffix_tree()
        sage:T.suffix_walk(((0,5),1))
        ("explicit", 0)
        sage:T.suffix_walk(((7,3),1))
        ("implicit", (9,4), 1)
    """
    #If the state is implicit
    parent=self.suffix_link(edge[0])
    for (i,j) in self._transition_function[edge[0]]:
        if self._transition_function[edge[0]][(i,j)]==edge[1]:
            break
    #(i-1,j) is the label of edge
    i-=1
    return self.count_and_skip(parent,(i,i+l))
    

from sage.combinat.words.suffix_trees import ImplicitSuffixTree
ImplicitSuffixTree.LZ_decomposition = LZ_decomposition
ImplicitSuffixTree.leftmost_covering_set = leftmost_covering_set
ImplicitSuffixTree.count_and_skip=count_and_skip
ImplicitSuffixTree.suffix_walk=suffix_walk

#===============================================================================
#        Fonction à ajouter Word
#===============================================================================
def longest_forward_extension(self,x,y):
    r"""
    Compute the length of le longest factor of self that starts at x and that
    matches a factor that starts at y. Returns 0 if x or y are not valid
    position in self.

    INPUTS:

        x,y - positions in self

    EXAMPLES:

        sage:w=Word('0011001')
        sage:w.longest_forward_extension(0,5)
        3
        sage:w.longest_forward_extension(0,2)
        0
        sage:w.longest_forward_extension(-3,2)
        0
    """
    length=self.length()
    if not (0<=x and 0<=y):
        return 0
    l=0
    while x<length and y<length and self[x]==self[y]:
        l+=1
        x+=1
        y+=1
    return l


def longest_backward_extension(self,x,y):
    r"""
    Compute the length of le longest factor of w that ends at x and that
    matches a factor that ends at y. Returns 0 if x or y are note valid position
    in self.

    INPUTS:

        x,y - positions in self

    EXAMPLES:
        
        sage:w=Word('0011001')
        sage:w.longest_backward_extension(7,2)
        3
        sage:w.longest_backward_extension(1,5)
        1
        sage:w.longest_forward_extension(4,23)
    """
    length=self.length()
    if not (x<length and y<length):
        return 0
    l=0
    while x>=0 and y>=0 and self[x]==self[y]:
        l+=1
        x-=1
        y-=1
    return l

from sage.combinat.words.finite_word import FiniteWord_class  
FiniteWord_class.longest_forward_extension=longest_forward_extension
FiniteWord_class.longest_backward_extension=longest_backward_extension

#===============================================================================
#        Nouvelle classe
#===============================================================================

class DecoratedSuffixTree(ImplicitSuffixTree):
    
    def __init__(self, w):
        r"""
        Construct the decorated suffix tree of a word

        A decorated suffix tree of w is the suffix tree of w  marked with the
        end point of all squares in the w.

        The symbol "$" is append to w to ensure de that each final state is a
        leaf of the suffix tree.

        When using pair as output, all the algorithm are linear in the length of
        the word w

        INPUT:

            w - word

        EXAMPLES:

            sage: w=Word('0011001')
            sage: DecoratedSuffixTree(w)
            Decorated suffix tree of : 0011001$

        REFERENCE:
        
          * [1] Gusfield, D., & Stoye, J. (2004). Linear time algorithms for 
          finding and representing all the tandem repeats in a string. Journal 
          of Computer and System Sciences, 69(4), 525-546.
        """
        if not isinstance(w, FiniteWord_class):
            raise ValueError("w must be a member of FiniteWord_class")
        if "$" in w:
            raise ValueError("The symbol '$' is reserved for this class ")
        end_symbol="$"
        w = Word(str(w)+"$")
        ImplicitSuffixTree.__init__(self, w)
        self.labeling=self._complete_labeling()

    def __repr__(self):
        w=self.word()
        if len(w)>40:
            w=str(w[:40])+'...'
        return "Decorated suffix tree of : %s" %w

    def _partial_labeling(self):
        r"""
        Make a depth first search in the suffix tree and mark some squares of a
        leftmost covering set of the tree. Used by _complete_labeling.

        EXAMPLES:

            sage:w=Word('abaababbabba')
            sage:T=DecoratedSuffixTree(w)
            sage:T._partial_labeling()
            {(3, 4): [1], (5, 1): [3], (5, 6): [1], (11, 17): [1], (13, 8): [1], (15, 10): [2]}
        """
        def node_processing(node,parent,(i,pos)):
            r"""
            Mark point along the edge (parent,node) if the string depht of parent is
            smaller than the lenght of the tamdem repeat at the head of P(node).
            Make it for all such squares pairs and remove them from P(node).
            P(node)=P[i][pos:]
            INPUTS:
                node - a node of T
                parent - the parent of node in T
                (i,pos) - the pair that represent the head of the list P(node)
            OUTPUT:
                (i,pos) - the new head of P(node)
            """
            while pos<len(P[i]) and P[i][pos]>string_depth[parent]:
                label=P[i][pos]-string_depth[parent]
                try:
                    labeling[(parent,node)].append(label)
                except KeyError:
                    labeling[(parent,node)]=[label]
                pos+=1
            return (i,pos)

        def treat_node(current_node,parent):
            r"""
            Proceed to a depth first search in T, couting the string_depth of
            each node a processing each node for marking

            To initiate de depth first search call treat_node(0,None)

            INPUTS:

                current_node - A node
                parent - Parent of current_node in T

            OUTPUT:

                The resulting list P(current_node) avec current_node have been
                process by node_processing. The ouput is a pair (i,pos) such
                that P[i][pos:] is the list of current_node.
            """
            #Call recursively on children of current_node
            if D.has_key(current_node):
                node_list=(n,0)
                for child in D[current_node].iterkeys():
                    (i,j)=D[current_node][child]
                    if j==None:
                        j=n
                    string_depth[child]=string_depth[current_node]+j-i
                    child_list=treat_node(child,current_node)
                    if child_list[0]<node_list[0]:
                        node_list=child_list
            else: #The node is a child
                node_list=(n-string_depth[current_node],0)
            #Make teatement on current node hear
            return node_processing(current_node,parent,node_list)
        
        P=leftmost_covering_set(self)
        D=self.transition_function_dictionary()
        string_depth=dict([(0,0)])
        n=len(self.word())
        labeling=dict()
        treat_node(0,None)
        return labeling

    def _complete_labeling(self):
        r"""
        Returns a dictionnary of edges of self, with markpoint for the end of
        each square types of T.word()

        INPUT:

            self - Suffix tree

        EXAMPLES:
            sage:w=Word('aabbaaba')
            sage:DecoratedSuffixTree(w)._complete_labeling()
            {(2, 7): [1], (5,4): [1]}
        """

        def walk_chain(u,v,l,start):
            r"""
            Execute a chain of suffix walk until a walk is unsuccesful or it got
            to a point already register in QP. Register all visited point in Q.

            INPUTS:

                (u,v) - edge on wich the point is registered
                l - depth of the registered point on (u,v)
                start - start of the squares registered by the label (u,v),l
            """
            #Mark the point in labeling
            try:
                labeling[(u,v)].append(l)
            except KeyError:
                labeling[(u,v)]=[l]
            #Make the walk
            final_state=self.suffix_walk(((u,v),l))
            successful=False
            if final_state[0]=='explicit':
                parent=final_state[1]
                transition=self._find_transition(parent,self._letters[start])
                if transition!=None:
                    child=transition[1]
                    successful=True
                    depth=1
            else:
                parent=final_state[1][0]
                child=final_state[1][1]
                depth=final_state[2]
                next_letter=self._letters[D[parent][child][0]+depth]
                if next_letter==self._letters[start]:
                    successful=True
                    depth+=1
            #If needed start a new walk
            if successful:
                try:
                    if depth not in prelabeling[(parent,child)]:
                        walk_chain(parent,child,depth,start+1)
                except KeyError:
                    walk_chain(parent,child,depth,start+1)

        def treat_node(current_node,(i,j)):
            r"""
            Execute a depht first search on self and start a suffix walk for
            labeled points on each edges of T. The fonction is reccursive, call
            treat_node(0,(0,0)) to initiate the search

            INPUTS:

                current_node - The node that is to treat
                (i,j) - Pair of index such that the path from 0 to current_node
                reads T.word()[i:j]
            """
            if D.has_key(current_node):
                for child in D[current_node].iterkeys():
                    edge=(current_node,child)
                    edge_label=D[edge[0]][edge[1]]
                    treat_node(child,(edge_label[0]-(j-i),edge_label[1]))
                    if prelabeling.has_key((current_node,child)):
                        for l in prelabeling[edge]:
                            square_start=edge_label[0]-(j-i)
                            walk_chain(current_node,child,l,square_start)

        prelabeling=self._partial_labeling()
        labeling=dict()
        D=self.transition_function_dictionary()
        treat_node(0,(0,0))
        return labeling
        
    def square_vocabulary(self,output="pair"):
        r"""
        Return the list of squares in the squares vocabulary of self.word.
        Return a list of pair in output="pair" and the explicit word if
        output="word"

        INPUTS:

            output - "pair" or "word"

        EXAMPLES:

            sage: w=Word('aabb')
            sage: DecoratedSuffixTree(w).square_vocabulary()
            [(0, 0), (0, 2), (2, 2)]
            sage: w=Word('00110011010')
            sage: DecoratedSuffixTree(w).square_vocabulary(output="word")
            [word: , word: 01100110, word: 00110011, word: 00, word: 11, word: 1010]
        """
        def treat_node(current_node,(i,j)):
            if D.has_key(current_node):
                for child in D[current_node].iterkeys():
                    edge=(current_node,child)
                    edge_label=(D[edge[0]][edge[1]])
                    treat_node(child,(edge_label[0]-(j-i),edge_label[1]))
                    if Q.has_key((current_node,child)):
                        for l in Q[(current_node,child)]:
                            square_start=edge_label[0]-(j-i)
                            pair=(square_start,edge_label[0]+l-square_start)
                            squares.append(pair)

        if not(output=="pair" or output=="word"):
            raise ValueError("output should be 'pair' or 'word'; got %s" %output)
        D=self.transition_function_dictionary()
        Q=self.labeling
        squares=[(0,0)]
        treat_node(0,(0,0))
        if output=="pair":
            return squares
        else:
            return [self.word()[i:i+l] for (i,l) in squares]

#===============================================================================
#        Pas encore triée
#===============================================================================

def naive_square_voc(T):
    r"""Compute the square vocabulary of T.word()
    INPUTS:
        T - Suffix Tree
    OUTPUT:
        Square vocabulary of T.word()"""
    squares=[]
    for f in [v for v in T.factor_iterator() if v.is_square()]:
        squares.append(f)
    return set(squares)

def run_test(n,alphabet='01',test_for_double=False):
    r"""
    Compare the algorithme with a naive algorithme
    INPUTS:
        n - size of words to test on
        test_for_double - If true detect a bug if the algorithm return a double
        alphabet - the alphabet to test on
    OUTPUT:
        True if works for all words, False if it bugs for a word
    """
    for w in tqdm(Words(alphabet,n)):
        S1=naive_square_voc(w)
        T=DecoratedSuffixTree(w)
        L=T.square_vocabulary(output="word")
        S2=set(L)
        if test_for_double and len(S2)!=len(L):
            print "Problème de doublon avec %s" %w
            print "is_sort_leftmost(%s)=%s" %(w,is_sort_leftmost(w))
            return False
        if not(S2.issubset(S1) and S1.issubset(S2)):
            print "Problème avec %s" %w
            return False
    return True

def is_sort_leftmost(w,strictly=True):
    r"""
    Verify if the leftmost covering set is sort according to the algorithm.
    If strictly is true, the list must be stricly increasing in order
    INPUTS:
        w - A word ending with $
    """
    covering=leftmost_covering_set(w.suffix_tree())
    for l in covering:
        for i in range(len(l)-1):
            if strictly and l[i][1]<=l[i+1][1]:
                return False
            elif l[i][1]<=l[i+1][1]:
                return False
    return True

def LZ_decomposition_explicit(T):
    r"""Take the explicit suffix tree of a word and return the Lempel-Ziv
    decomposition of the word in the form of a list iB of index such that the
    blocks of the decomposition are T.word()[iB[k]:iB[k+1]]"""
    iB=[0]
    i=0
    w=T.word()
    while i<len(T.word()):
        l=0
        s=0
        ((x,y),successor)=T._find_transition(s, w[i])
        x=x-1 #Pourquoi find_transtion retourne pas la bonne étiquette?
        while x<i+l and y!=None and i+l<len(w):
            l+=y-x
            s=successor
            if i+l==len(w):
                break
            transition=T._find_transition(s,w[i+l])
            ((x,y),successor)=transition
            x=x-1 #ici
        i+=max(1,l)
        iB.append(i)
    return iB
