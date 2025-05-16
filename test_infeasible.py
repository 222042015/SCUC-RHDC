import gurobipy as gp

model = gp.read("rh_cb.lp")

model.computeIIS()

# gp.write("rh_cb.ilp", model.IIS)
model.write("rh_cb.ilp")