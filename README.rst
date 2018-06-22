Linear algorithm to find square vocabulary of words
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This repository contains code for computing all the different square factors of
a words. This program is an implementation of the algorithm described of the
article *Linear time algorithms for finding and representing all the tandem
repeats in a string* by Dan Gusfield and Jens Stoye
(https://doi.org/10.1016/j.jcss.2004.03.004)

A detailed report (in french) is available in the repo.

Background
==========

Given a word `u` we say that `u` is a square if `u=vv`. The square vocabulary of
a words `w=w_0w_1...w_n` is the set of all the distinct square that occur as a
factor in `w`. Since a factor of length `l` starting at position `i` can be
represent by a pair `(i,l)`, the square vocabulary is a set of such pair.
It has been shown that the number of factors in the square vocabulary is in
`O(n)`.

The algorithm is based on the `Suffix tree
<https://en.wikipedia.org/wiki/Suffix_tree>`__ data structure.

Dependencies
============

One must have `Sagemath <http://www.sagemath.org>`__ installed to
run the program. It is mostly used to for the words and the suffix tree
library.

How to use
==========

The main class is ``DecoratedSuffixTree``. A instance of this class is used to
compute compute the square vocabulary of a words in the form of a list of pair
`(i,l)`.

Below are some examples that can be reproduced once Sagemath is started ::

    sage: load('decorated_suffix_tree.py')
    sage: w=Word('aabb')
    sage: DecoratedSuffixTree(w).square_vocabulary()
    [(0, 0), (0, 2), (2, 2)]
    sage: w=Word('00110011010')
    sage: DecoratedSuffixTree(w).square_vocabulary(output="word")
    [word: , word: 01100110, word: 00110011, word: 00, word: 11, word: 1010]

License
=======

All files in this repository are subject to the `GPLv3 license
<https://www.gnu.org/licenses/gpl-3.0.en.html>`__.
