
class Test():
    
    def __init__(self, left, right):
        self.left = left
        self.right = right
        
    def update_right(self):
        self.right.h += 1
        self.right.u += 1
        
    def update_left(self):
        self.left.h += 2
        self.left.u += 1
        

class obj():
    
    def __init__(self, h, u):
        self.h = h
        self.u = u
        
    def reset_h(self):
        self.h = 0
        
    def reset_u(self):
        self.u = 0
        
left = obj(1,1)
mid = obj(0,0)
right = obj(1,1)
cells = [Test(left, mid), Test(mid, right)]

cells[0].update_right()
print(cells[0].right.h, cells[1].left.h)
cells[1].update_left()
print(cells[0].right.h, cells[1].left.h)
mid.reset_h()
print(cells[0].right.h, cells[1].left.h)