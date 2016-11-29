# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: test.py
@time: 11/28/2016 7:06 PM
"""

from learningPY.test1 import my_abs, f1, normalize

a = my_abs(11)
print(a)

#print(my_abs("d"))

b = (1,2,3,4,5,6,7)
kw = {'A':11,'B':12}
print(f1(*b,**kw))

print(list(map(str, [1, 2, 3, 4, 5, 6, 7, 8, 9])))
print({'1':'a','2':'b','3':'c'}['2'])

#map  &  reduce
L1 = ['adam', 'LISA', 'barT', 'jiNCHeNg']
L2 = list(map(normalize, L1))
print(L2)