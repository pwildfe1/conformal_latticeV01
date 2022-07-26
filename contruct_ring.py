import rhinoscriptsyntax as rs
import math as m

def main():
    section = rs.GetObject("please select section", rs.filter.curve)
    sections = []
    rs.EnableRedraw(False)
    for i in range(500):
        sections.append(rs.RotateObject(section, [0, 0, 0], i/500 * 360, copy = True))
    grid = []
    for i in range(len(sections)):
        divPts = rs.DivideCurve(sections[i], 200)
        for j in range(len(divPts)):
            grid.append(divPts[j])
    srf = rs.AddSrfPtGrid([500, 200], grid, closed = (True, True))

main()