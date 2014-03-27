"""
http://rhettinger.wordpress.com/2011/05/26/super-considered-super/
"""

class Combination(list):
	def __init__(self,n,m,a,b):
		super(Combination,self).__init__()
		self._n,self._m = n,m
		self._first,self._last = a,b
		self._count = 0
		self._generate()

	def _add(self,l):
		if self._count >= self._first and self._count <= self._last:
			self.append(l)
		self._count += 1

	def _generate(self,b=1,l=[0]):
		n,m = self._n,self._m
		L = len(l)
		if L == m+1:
			self._add(l[1:])
		else:
			c = l[-1]
			q,r = c+1,min(c+b+m,n+1)
			for i in range(q,r):
				self._generate(b+1,l+[i])
