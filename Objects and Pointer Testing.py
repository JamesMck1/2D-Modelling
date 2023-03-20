
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
print(cells[0].right.h, cells[1].left.h, mid.h)
cells[1].update_left()
print(cells[0].right.h, cells[1].left.h, mid.h)
mid.reset_h()
print(cells[0].right.h, cells[1].left.h, mid.h)

class Cell():
    
    def __init__(self, h, u):
        self.h = h
        self.u = u
        
    def assign_interface(self, int_L, int_R):
        self.int_L = int_L
        self.int_R = int_R
        
    def update(self):
        self.h = self.h + 0.5*(self.int_L.f1 - self.int_R.f1)
        self.u = self.h + 0.5*(self.int_L.f2 - self.int_R.f2)

class Clones():
    
    def __init__(self, Cell):
        self.Cell = Cell
        self.h = self.Cell.h
        self.u = self.Cell.u
    
    #def __getattr__(self, attr):
    #    return getattr(self.Cell, attr)
    
    def update(self):
        self.h = self.Cell.h
        self.u = self.Cell.u
    
class Clones_2():
    
    def __init__(self, Cell):
        self.h = Cell.h
        self.u = Cell.u

class Interface():
    
    def __init__(self, left_cell, right_cell):
        self.left_cell = left_cell
        self.right_cell = right_cell
        
    def fluxes(self):
        self.f1 = self.left_cell.h - self.right_cell.h
        self.f2 = self.left_cell.u - self.right_cell.u
        
cell_1 = Cell(4,0)
cell_2 = Cell(4,2)
cell_3 = Cell(1,3)
Clone = Clones(cell_2)
Clone_2 = Clones_2(cell_2)
Int_1 = Interface(cell_1, cell_2)
Int_2 = Interface(cell_2, cell_3)
cell_2.assign_interface(Int_1, Int_2)
Int_1.fluxes()
Int_2.fluxes()
cell_2.update()

print(Clone.h, Clone.u)

Clone.update()

print(Clone.Cell.h, Clone_2.h, Clone.h)    
    