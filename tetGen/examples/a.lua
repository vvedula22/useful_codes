-- Lua script.
p=tetview:new()
p:load_plc("brain.smesh")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
