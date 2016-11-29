# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: test1.py
@time: 11/28/2016 7:11 PM
"""
"""函数的参数：
1.位置参数
2.默认参数：指向不变对象；
3.可变参数，i.用list/tuple；ii.*[1,2,3]
4.关键字参数：
>>> extra = {'city': 'Beijing', 'job': 'Engineer'}
>>> person('Jack', 24, **extra)
name: Jack age: 24 other: {'city': 'Beijing', 'job': 'Engineer'}
def person(name, age, **kw):
    print('name:', name, 'age:', age, 'other:', kw)
**代表吧dist以key-value的形式传入

对于任意函数，都可以通过类似func(*args, **kw)的形式调用它，
无论它的参数是如何定义的。
"""
def my_abs(x=11):   #默认参数必须指向不变对象！
    if not isinstance(x,(int,float)):
        raise TypeError('bad operand type')
    if x > 0:
        return x
    if x <= 0:
        return -x

def nop():
    pass

def f1(a, b, c=11, *d, **kw):
    return 'a=',a,'b=',b,'c=', c,'d=', d,'kw=', kw

'''切片List[:10:2]前十个数，每两个取一个
列表生成器 [a*a for a in range(1,100)]
generator ()'''


'''map/reduce:
'''

def normalize(name):
    return name[:1].upper()+name[1:].lower()