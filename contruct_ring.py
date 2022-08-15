import rhinoscriptsyntax as rs
import math as m

def main():
    section = rs.GetObject("please select section", rs.filter.curve)
    sections = []
    rs.EnableRedraw(False)
    resoU, resoV = 300, 50
    for i in range(resoU):
        sections.append(rs.RotateObject(section, [0, 0, 0], i/(resoU-1) * 360, copy = True))
    grid = []
    for i in range(len(sections)):
        divPts = rs.DivideCurve(sections[i], resoV)
        for j in range(len(divPts)):
            grid.append(divPts[j])
    srf = rs.AddSrfPtGrid([resoU, resoV], grid, closed = (True, True))

main()