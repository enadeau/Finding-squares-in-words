from tqdm import tqdm
#To install tqdm run sage -pip install tqdm in a terminal

def LZ_decomposition(T):
    r"""Take the implicit suffix tree of a word and return the Lempel-Ziv
    decomposition of the word in the form of a list iB of index such that the
    blocks of the decomposition are T.word()[iB[k]:iB[k+1]]"""
    iB=[0]
    i=0
    w=T.word()
    while i<len(w):
        l=0
        s=0
        ((x,y),successor)=T._find_transition(s, w[i])
        x=x-1 #Pourquoi find_transtion retourne pas la bonne étiquette?
        while x<i+l:
            if y==None:
                l=len(w)-i
            else:
                l+=y-x
            if i+l==len(w):
                break
            s=successor
            ((x,y),successor)=T._find_transition(s,w[i+l])
            x=x-1 #même question ici
        i+=max(1,l)
        iB.append(i)
    return iB

def leftmost_covering_set(T):
    r"""Compute the leftmost covering set of squares in T.word().
    INPUTS:
        T - Suffix tree
    OUTPUTS:
        Leftmost covering set of pair
    """
    def condition1_square_pairs(i):
        r"""Compute the square that has their center in the i-th block of 
        LZ-decomposition and that start in the i-th block and end in the
        (i+1)-th"""
        for k in range(1,B[i+1]-B[i]+1):
            q=B[i+1]-k
            k1=longest_forward_extension(w,B[i+1],q)
            k2=longest_backward_extension(w,B[i+1]-1,q-1)
            start=max(q-k2,q-k+1)
            if k1+k2>=k and k1>0:# and start>=B[i]:
                #print "Condition 1 yield (%s,%s) for block %s" %(start,2*k,i)
                #print start>=B[i]
                yield (start,2*k)
    
    def condition2_square_pairs(i):
        r"""Compute the squares that has their center in the i-th block of the
        LZ-decomposition and that starts in the (i-1)-th block or before. Their
        end is either in the i-th or the (i+1)-th block"""
        try:
            end=B[i+2]-B[i]+1
        except IndexError:
            end=B[i+1]-B[i]+1
        for k in range(1,end):
            q=B[i]+k
            k1=longest_forward_extension(w,B[i],q)
            k2=longest_backward_extension(w,B[i]-1,q-1)
            start=max(B[i]-k2,B[i]-k+1)
            #print "k=%s q=%s k1=%s k2=%s, start=%s, h1=%s" %(k,q,k1,k2,start,B[i+1])
            if k1+k2>=k and k1>0 and start+k<=B[i+1] and k2>0:# and start<B[i]:
                #print "Condition 2 yield (%s,%s) for block %s" %(start,2*k,i)
                #print start<B[i]
                yield (start,2*k)

    w=T.word()
    B=LZ_decomposition(T)
    P=[[] for _ in w]
    for i in range(len(B)-1):
        squares=list(condition2_square_pairs(i))+list(condition1_square_pairs(i))
        for (i,l) in squares:
            P[i].append((i,l))
    for l in P:
        l.reverse()
    return P

def partial_labeling(T):
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
        while pos<len(P[i]) and P[i][pos][1]>string_depth[parent]:
            label=P[i][pos][1]-string_depth[parent]
            try:
                labeling[(parent,node)].append(label)
            except KeyError:
                labeling[(parent,node)]=[label]
            pos+=1
        return (i,pos)

    def treat_node(current_node,parent):
        r"""
        Proceed to a depth first search in T, couting the string_depth of each
        node a processing each node for marking

        To initiate de depth first search call treat_node(0,None)

        INPUTS:
            current_node - A node
            parent - Parent of current_node in T
        OUTPUT:
            The resultint list P(current_node) avec current_node have been
            process by node_processing. The ouput is a pair (i,pos) such that
            P[i][pos:] is the list of current_node
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
    
    P=leftmost_covering_set(T)
    D=T.transition_function_dictionary()
    string_depth=dict([(0,0)])
    n=len(T.word())
    labeling=dict()
    treat_node(0,None)
    return labeling

def complete_labeling(T):
    r"""
    Returns a dictionnary of edges of T, whit markpoint for the end of each
    square types of T.word()
    INPUT:
        T - Suffix tree
    """

    def count_and_skip(node,(i,j)):
        r"""
        Use count and skip trick to follow the path starting and "node" and
        reading T.word()[i:j]. We assume that reading T.word()[i:j] is possible
        from "node"

        INPUTS:
            node - explicit node of T
            (i,j) - Indices of factor T.word()[i:j]
        OUTPUT:
            The node obtained by starting at "node" and following the edges
            labeled by the letter of T.word()[i:j]. Returns "("explicit",
            end_node) if w and at a "end_node", of "("implicit", (edge, d))" if
            we end at d sports along an edge.
        """
        if i==j: #We're done reading the factor
            return ('explicit',node)
        transition=T._find_transition(node,T._letters[i])
        child=transition[1]
        if transition[0][1]==None: #The child is a leaf
            edge_length=len(T.word())-transition[0][0]+1
        else:
            edge_length=transition[0][1]-transition[0][0]+1
        if edge_length>j-i: #The reading stop on this edge
            return ('implicit',((node,child),j-i))
        return count_and_skip(child,(i+edge_length,j))

    def suffix_link_walk(u,v,l,start):
        r"""
        Execute a chain of suffix walk until a walk is unsuccesful or it got to
        a point already register in QP. Register all visited point in Q
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
        (i,j)=(D[u][v][0],D[u][v][0]+l) #suffix to read
        u=T.suffix_link(u)
        final_state=count_and_skip(u,(i,j))
        successful=False
        if final_state[0]=='explicit':
            parent=final_state[1]
            transition=T._find_transition(parent,T._letters[start])
            if transition!=None:
                child=transition[1]
                successful=True
                depth=1
        else:
            parent=final_state[1][0][0]
            child=final_state[1][0][1]
            depth=final_state[1][1]
            next_letter=T._letters[D[parent][child][0]+depth]
            if next_letter==T._letters[start]:
                successful=True
                depth+=1
        #If needed start a new walk
        if successful:
            try:
                if depth not in prelabeling[(parent,child)]:
                    suffix_link_walk(parent,child,depth,start+1)
            except KeyError:
                suffix_link_walk(parent,child,depth,start+1)

    def treat_node(current_node,(i,j)):
        r"""
        Execute a depht first search on T and start a suffix walk for labeled
        points on each edges of T. The fonction is reccursive, call
        treat_node(0,(0,0)) to start the search
        INPUTS
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
                        suffix_link_walk(current_node,child,l,square_start)

    prelabeling=partial_labeling(T)
    labeling=dict()
    D=T.transition_function_dictionary()
    treat_node(0,(0,0))
    return labeling 

def list_squares(T):
    def treat_node(current_node,(i,j)):
        if D.has_key(current_node):
            for child in D[current_node].iterkeys():
                edge=(current_node,child)
                edge_label=(D[edge[0]][edge[1]])
                treat_node(child,(edge_label[0]-(j-i),edge_label[1]))
                if Q.has_key((current_node,child)):
                    for l in Q[(current_node,child)]:
                        square_start=edge_label[0]-(j-i)
                        squares.append(T.word()[square_start:edge_label[0]+l])
    D=T.transition_function_dictionary()
    Q=complete_labeling(T)
    squares=[Word('')]
    treat_node(0,(0,0))
    return squares

#Truc à améliorer
def longest_forward_extension(w,x,y):
    r"""Compute the length of le longest factor of w that start at x and that
    matches a factor that start at y.
    INPUTS:
        w - a word
        x,y - position in w
    OUTPUTS:
        Length of the longest foward extension"""
    l=0
    while x<len(w) and y<len(w) and w[x]==w[y]:
        l+=1
        x+=1
        y+=1
    return l


def longest_backward_extension(w,x,y):
    r"""Compute the length of le longest factor of w that ends at x and that
    matches a factor that ends at y.
    INPUTS:
        w - a word
        x,y - position in w
    OUTPUTS:
        Length of the longest backward extension"""
    l=0
    while x>=0 and y>=0 and w[x]==w[y]:
        l+=1
        x-=1
        y-=1
    return l

#Truc qui serve à rien
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
        w=Word(w,alphabet+'$')*Word('$')
        T=w.suffix_tree()
        L=list_squares(T)
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

def condition1squares(T):
    B=LZ_decomposition(T)
    for i in range(len(B)-1):
        for k in range(1,B[i+1]-B[i]+1):
            q=B[i+1]-k
            k1=longest_forward_extension(T.word(),B[i+1],q)
            k2=longest_backward_extension(T.word(),B[i+1]-1,q-1)
            start=max(q-k2,q-k+1)
            if k1+k2>=k and k1>0:
                yield (start,2*k)

def condition2squares(T):
    B=LZ_decomposition(T)
    for i in range(len(B)-1):
        try:
            end=B[i+2]-B[i]+1
        except IndexError:
            end=B[i+1]-B[i]+1
        for k in range(1,end):
            q=B[i]+k
            k1=longest_forward_extension(T.word(),B[i],q)
            k2=longest_backward_extension(T.word(),B[i]-1,q-1)
            start=max(B[i]-k2,B[i]-k+1)
            if k1+k2>=k and k1>0 and start+k<=B[i+1] and k2>0:
                yield (start,2*k)
