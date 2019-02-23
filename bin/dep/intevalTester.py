from interval import Interval, IntervalSet
r1 = IntervalSet([Interval(1,1000),Interval(1001,2000),Interval(2001,3000)])
r2 = IntervalSet([Interval(50,52),Interval(700,750),Interval(1125,1444)])


#r1 = IntervalSet([Interval(1, 1000), Interval(1100, 1200)])
#r2 = IntervalSet([Interval(30, 50), Interval(60, 200), Interval(1150, 1300)])

for k in  r1 - r2:
	print k
