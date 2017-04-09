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
            if k1+k2>=k and k1>0 and start>=B[i]:
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
            if k1+k2>=k and k1>0 and start+k<=B[i+1] and k2>0:
                yield (start,2*k)

    w=T.word()
    B=LZ_decomposition(T)
    P=[[] for _ in w]
    for i in range(len(B)-1):
        squares=list(condition1_square_pairs(i))+list(condition2_square_pairs(i))
        for (i,l) in squares:
            P[i].append((i,l))
    for l in P:
        l.reverse()
    return P

def partial_labeling(T):
    def treat_node(current_node,parent):
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
        (i,pos)=node_list
        while len(P[i])>pos and P[i][pos][1]>string_depth[parent]:
            if edge_label.has_key((parent,current_node)):
                edge_label[(parent,current_node)].append(P[i][pos][1]-string_depth[parent])
            else:
                edge_label[(parent,current_node)]=[P[i][pos][1]-string_depth[parent]]
            pos+=1
        return (i,pos)
    
    P=leftmost_covering_set(T)
    D=T.transition_function_dictionary()
    string_depth=dict([(0,0)])
    n=len(T.word())
    edge_label=dict()
    treat_node(0,None)
    return edge_label

def complete_labelling(T):
    def count_and_skip(T,node,(i,j)):
        r"""
        Use count and skip trick to follow the path starting and "node" and
        reading T.word()[i:j]. We assume that reading T.word()[i:j] is possible
        from "node"

        INPUTS:
            T - suffix tree of a word w
            node - explicit node of T
            (i,j) - Indices of factor T.word()[i:j]
        OUTPUT:
            The node obtained by starting at "node" and following the edges
            labelled by the letter of T.word()[i:j]. Returns "("explicit",
            end_node) if w and at a "end_node", of "("implicit", (edge, d))" if
            we end at d sports along an edge. If it's not possible to read the
            factor, return None
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
        return count_and_skip(T,child,(i+edge_length,j))

    def suffix_link_walk(u,v,l,start):
        print (u,v,l,start)
        (i,j)=(D[u][v][0],D[u][v][0]+l) #suffix to read
        u=T.suffix_link(u)
        v=T._find_transition(u,T.word()[i])[1]
        while D[u][v][1]!=None and j-i>=D[u][v][1]-D[u][v][0]:
            l_edge=(D[u][v][1]-D[u][v][0])
            i=i+l_edge
            u=v
            v=T._find_transition(v,T.word()[i])[1]
        print (u,v,l,start)
        print (i,j)
        transition=T.transition_function(T.word()[i:j]*T.word()[start:start+1],node=u)
        if transition!=None:
            v=T._find_transition(u,(T.word()[i:j]*T.word()[start:start+1])[0])[1]
            has_transition=True
        else:
            has_transition=False
        is_not_registered=((not QP.has_key((u,v))) or (QP.has_key((u,v)) and (j-i+1) not in QP[(u,v)]))
        if has_transition and is_not_registered:
            if Q.has_key((u,v)):
                Q[(u,v)].append(j-i+1)
            else:
                Q[(u,v)]=[j-i+1]
            if (not QP.has_key((u,v))) or (QP.has_key((u,v)) and (j-i+1) not in QP[(u,v)]):
                suffix_link_walk(u,v,j-i+1,start+1)

    def treat_node(current_node,(i,j)):
        if D.has_key(current_node):
            for child in D[current_node].iterkeys():
                edge=(current_node,child)
                edge_label=(D[edge[0]][edge[1]])
                treat_node(child,(edge_label[0]-(j-i),edge_label[1]))
                if QP.has_key((current_node,child)):
                    for l in QP[(current_node,child)]:
                        square_start=edge_label[0]-(j-i)
                        suffix_link_walk(current_node,child,l,square_start)
    QP=partial_labeling(T)
    Q=dict()
    for key in QP.iterkeys():
        Q[key]=copy(QP[key])
    D=T.transition_function_dictionary()
    treat_node(0,(0,0))
    return Q

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
    Q=complete_labelling(T)
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

def run_test(n,print_word=False):
    for w in Words('012',n):
        if print_word==True:
            print w
        S1=naive_square_voc(w)
        w=Word(str(w)+'$')
        T=w.suffix_tree()
        L=list_squares(T)
        S2=set(L)
        if len(S2)!=len(L):
            return False
        if not(S2.issubset(S1) and S1.issubset(S2)):
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
