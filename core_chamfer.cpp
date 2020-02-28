#include <TopoDS_Shape.hxx> 
#include <TopoDS.hxx> 
#include <BRepPrimAPI_MakeBox.hxx> 
#include <TopoDS_Solid.hxx> 
#include <BRepFilletAPI_MakeFillet.hxx> 
#include <TopExp_Explorer.hxx> 

TopoDS_Shape FilletedBox(const Standard_Real a, 
                        const Standard_Real  b, 
                        const Standard_Real  c, 
                        const Standard_Real  r) 
{ 
    TopoDS_Solid Box =  BRepPrimAPI_MakeBox(a,b,c); 
    BRepFilletAPI_MakeFillet  MF(Box); 
    
    // add all the edges  to fillet 
    TopExp_Explorer  ex(Box,TopAbs_EDGE); 
    while (ex.More()) 
    { 
    MF.Add(r,TopoDS::Edge(ex.Current())); 
    ex.Next(); 
    } 
    return MF.Shape(); 
} 

void CSampleTopologicalOperationsDoc::OnEvolvedblend1() 
{ 
    TopoDS_Shape theBox  = BRepPrimAPI_MakeBox(200,200,200); 
    BRepFilletAPI_MakeFillet  Rake(theBox); 
    ChFi3d_FilletShape  FSh = ChFi3d_Rational; 
    Rake.SetFilletShape(FSh); 
    TColgp_Array1OfPnt2d  ParAndRad(1, 6); 
    ParAndRad(1).SetCoord(0.,  10.); 
    ParAndRad(1).SetCoord(50.,  20.); 
    ParAndRad(1).SetCoord(70.,  20.); 
    ParAndRad(1).SetCoord(130.,  60.); 
    ParAndRad(1).SetCoord(160.,  30.); 
    ParAndRad(1).SetCoord(200.,  20.); 
    TopExp_Explorer  ex(theBox,TopAbs_EDGE); 
    Rake.Add(ParAndRad, TopoDS::Edge(ex.Current())); 
    TopoDS_Shape  evolvedBox = Rake.Shape(); 
} 
