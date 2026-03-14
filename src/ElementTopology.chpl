/**
    This module contains the api to CGNS constants for elements type definition and
    procedures to convert the element type into number of facets and number of nodes.
**/
module ElementTopology
{
use globalParam;

/* Element types from CGNS standard. */
/** Null element type from cgnslib.h: Value=0 **/
const ElementTypeNull : c_int = 0;
/** User-defined element type from cgnslib.h: Value=1 **/
const ElementTypeUserDefined : c_int = 1;
/** 0-D node element type from cgnslib.h: Value=2 **/
const NODE : c_int = 2;
/** 1-D linear bar element type from cgnslib.h: Value=3 **/
const BAR_2 : c_int = 3;
/** 1-D quadratic bar element type from cgnslib.h: Value=4 **/
const BAR_3 : c_int = 4;
/** 2-D linear triangle element type from cgnslib.h: Value=5 **/
const TRI_3 : c_int = 5;
/** 2-D quadratic triangle element type from cgnslib.h: Value=6 **/
const TRI_6 : c_int = 6;
/** 2-D linear quadrangle element type from cgnslib.h: Value=7 **/
const QUAD_4 : c_int = 7;
/** 2-D quadratic quadrangle element type from cgnslib.h: Value=8 **/
const QUAD_8 : c_int = 8;
/** 2-D quadratic quadrangle element type from cgnslib.h: Value=9 **/
const QUAD_9 : c_int = 9;
/** 3-D linear tetrahedron element type from cgnslib.h: Value=10 **/
const TETRA_4 : c_int = 10;
/** 3-D quadratic tetrahedron element type from cgnslib.h: Value=11 **/
const TETRA_10 : c_int = 11;
/** 3-D linear pyramid element type from cgnslib.h: Value=12 **/
const PYRA_5 : c_int = 12;
/** 3-D quadratic pyramid element type from cgnslib.h: Value=13 **/
const PYRA_14 : c_int = 13;
/** 3-D linear pentahedron element type from cgnslib.h: Value=14 **/
const PENTA_6 : c_int = 14;
/** 3-D quadratic pentahedron element type from cgnslib.h: Value=15 **/
const PENTA_15 : c_int = 15;
/** 3-D quadratic pentahedron element type from cgnslib.h: Value=16 **/
const PENTA_18 : c_int = 16;
/** 3-D linear hexahedron element type from cgnslib.h: Value=17 **/
const HEXA_8 : c_int = 17;
/** 3-D quadratic hexahedron element type from cgnslib.h: Value=18 **/
const HEXA_20 : c_int = 18;
/** 3-D quadratic hexahedron element type from cgnslib.h: Value=19 **/
const HEXA_27 : c_int = 19;
/** Mixed element type from cgnslib.h: Value=20 **/
const MIXED : c_int = 20;
/** 3-D quadratic pyramid element type from cgnslib.h: Value=21 **/
const PYRA_13 : c_int = 21;
/** N-polygon element type from cgnslib.h: Value=22 **/
const NGON_n : c_int = 22;
/** N-polyhedron element type from cgnslib.h: Value=23 **/
const NFACE_n : c_int = 23;
/** 1-D cubic bar element type from cgnslib.h: Value=24 **/
const BAR_4 : c_int = 24;
/** 2-D cubic triangle element type from cgnslib.h: Value=25 **/
const TRI_9 : c_int = 25;
/** 2-D cubic triangle element type from cgnslib.h: Value=26 **/
const TRI_10 : c_int = 26;
/** 2-D cubic quadrangle element type from cgnslib.h: Value=27 **/
const QUAD_12 : c_int = 27;
/** 2-D cubic quadrangle element type from cgnslib.h: Value=28 **/
const QUAD_16 : c_int = 28;
/** 3-D cubic tetrahedron element type from cgnslib.h: Value=29 **/
const TETRA_16 : c_int = 29;
/** 3-D cubic tetrahedron element type from cgnslib.h: Value=30 **/
const TETRA_20 : c_int = 30;
/** 3-D cubic pyramid element type from cgnslib.h: Value=31 **/
const PYRA_21 : c_int = 31;
/** 3-D cubic pyramid element type from cgnslib.h: Value=32 **/
const PYRA_29 : c_int = 32;
/** 3-D cubic pyramid element type from cgnslib.h: Value=33 **/
const PYRA_30 : c_int = 33;
/** 3-D cubic pentahedron element type from cgnslib.h: Value=34 **/
const PENTA_24 : c_int = 34;
/** 3-D cubic pentahedron element type from cgnslib.h: Value=35 **/
const PENTA_38 : c_int = 35;
/** 3-D cubic pentahedron element type from cgnslib.h: Value=36 **/
const PENTA_40 : c_int = 36;
/** 3-D cubic hexahedron element type from cgnslib.h: Value=37 **/
const HEXA_32 : c_int = 37;
/** 3-D cubic hexahedron element type from cgnslib.h: Value=38 **/
const HEXA_56 : c_int = 38;
/** 3-D cubic hexahedron element type from cgnslib.h: Value=39 **/
const HEXA_64 : c_int = 39;
/** 1-D quartic bar element type from cgnslib.h: Value=40 **/
const BAR_5 : c_int = 40;
/** 2-D quartic triangle element type from cgnslib.h: Value=41 **/
const TRI_12 : c_int = 41;
/** 2-D quartic triangle element type from cgnslib.h: Value=42 **/
const TRI_15 : c_int = 42;
/** 2-D quartic quadrangle element type from cgnslib.h: Value=43 **/
const QUAD_P4_16 : c_int = 43;
/** 2-D quartic quadrangle element type from cgnslib.h: Value=44 **/
const QUAD_25 : c_int = 44;
/** 3-D quartic tetrahedron element type from cgnslib.h: Value=45 **/
const TETRA_22 : c_int = 45;
/** 3-D quartic tetrahedron element type from cgnslib.h: Value=46 **/
const TETRA_34 : c_int = 46;
/** 3-D quartic tetrahedron element type from cgnslib.h: Value=47 **/
const TETRA_35 : c_int = 47;
/** 3-D quartic pyramid element type from cgnslib.h: Value=48 **/
const PYRA_P4_29 : c_int = 48;
/** 3-D quartic pyramid element type from cgnslib.h: Value=49 **/
const PYRA_50 : c_int = 49;
/** 3-D quartic pyramid element type from cgnslib.h: Value=50 **/
const PYRA_55 : c_int = 50;
/** 3-D quartic pentahedron element type from cgnslib.h: Value=51 **/
const PENTA_33 : c_int = 51;
/** 3-D quartic pentahedron element type from cgnslib.h: Value=52 **/
const PENTA_66 : c_int = 52;
/** 3-D quartic pentahedron element type from cgnslib.h: Value=53 **/
const PENTA_75 : c_int = 53;
/** 3-D quartic hexahedron element type from cgnslib.h: Value=54 **/
const HEXA_44 : c_int = 54;
/** 3-D quartic hexahedron element type from cgnslib.h: Value=55 **/
const HEXA_98 : c_int = 55;
/** 3-D quartic hexahedron element type from cgnslib.h: Value=56 **/
const HEXA_125 : c_int = 56;

/** Domain for arrays using element types as index **/
const elemTypeRange : range = ElementTypeNull..HEXA_125;

/** Array defining the number of nodes for each element type **/
const numNodesForElementType : [elemTypeRange] int = [-1, 0, 1, 2, 3, 3, 6, 4, 8, 9, 4, 10, 5, 14, 6, 15, 18, 8, 20, 27, 0, 13, 0, 0, 4, 9, 10, 12, 16, 16, 20, 21, 29, 30, 24, 38, 40, 32, 56, 64, 5, 12, 15, 16, 25, 22, 34, 35, 29, 50, 55, 33, 66, 75, 44, 98, 125];

/** Basic 0-D element defining a vertex. **/
proc Node_t()
{
    return (NODE, "NODE", 1, 0, 0, [ElementTypeNull], [-1], [-1], [ElementTypeNull], [-1], [-1]);
}

/** Basic 1-D element defining a linear edge. **/
proc Bar2_t()
{
    return (BAR_2, "BAR_2", 2, 1, 0,
                    [BAR_2],
                    [0, 2],
                    [0, 1],
                    [ElementTypeNull],
                    [-1],
                    [-1]);
}

/** 2-D element defining a linear triangular facet. **/
proc Tri3_t()
{
    return (TRI_3, "TRI_3", 3, 3, 1,
                    [BAR_2, BAR_2, BAR_2],
                    [0, 2, 4, 6],
                    [0, 1, 1, 2, 2, 0],
                    [TRI_3],
                    [0, 3],
                    [0, 1, 2]);
}

/** 2-D element defining a linear quadrangular facet. **/
proc Quad4_t()
{
    return (QUAD_4, "QUAD_4", 4, 4, 1,
                    [BAR_2, BAR_2, BAR_2, BAR_2],
                    [0, 2, 4, 6, 8],
                    [0, 1, 1, 2, 2, 3, 3, 0],
                    [QUAD_4],
                    [0, 4],
                    [0, 1, 2, 3]);
}

/** 3-D element defining a linear tetrahedral cell. **/
proc Tetra4_t()
{
    return (TETRA_4, "TETRA_4", 4, 6, 4,
                    [BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2],
                    [0, 2, 4, 6, 8, 10, 12],
                    [0, 1, 1, 2, 2, 0, 0, 3, 1, 3, 2, 3],
                    [TRI_3, TRI_3, TRI_3, TRI_3],
                    [0, 3, 6, 9, 12],
                    [0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3]);
 }

/** 3-D element defining a linear pyramidal cell. **/
proc Pyra5_t()
{
    return (PYRA_5, "PYRA_5", 5, 8, 5,
                    [BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2],
                    [0, 2, 4, 6, 8, 10, 12, 14, 16],
                    [0, 1, 1, 2, 2, 3, 3, 0, 0, 4, 1, 4, 2, 4, 3, 4],
                    [QUAD_4, TRI_3, TRI_3, TRI_3, TRI_3],
                    [0, 4, 7, 10, 13, 16],
                    [0, 3, 2, 1, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4]);
 }

/** 3-D element defining a linear prism cell. **/
proc Penta6_t()
{
    return (PENTA_6, "PENTA_6", 6, 9, 5,
                    [BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2],
                    [0, 2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [0, 1, 1, 2, 2, 0, 0, 3, 1, 4, 2, 5, 3, 4, 4, 5, 5, 3],
                    [QUAD_4, QUAD_4, QUAD_4, TRI_3, TRI_3],
                    [0, 4, 8, 12, 15, 18],
                    [0, 1, 4, 3, 1, 2, 5, 4, 2, 0, 3, 5, 0, 2, 1, 3, 4, 5]);
}

/** 3-D element defining a linear hexahedral cell. **/
proc Hexa8_t()
{
    return (HEXA_8, "HEXA_8", 8, 12, 6,
                    [BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2, BAR_2],
                    [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24],
                    [0, 1, 1, 2, 2, 3, 3, 0, 0, 4, 1, 5, 2, 6, 3, 7, 4, 5, 5, 6, 6, 7, 7, 4],
                    [QUAD_4, QUAD_4, QUAD_4, QUAD_4, QUAD_4, QUAD_4],
                    [0, 4, 8, 12, 16, 20, 24],
                    [0, 3, 2, 1, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 0, 4, 7, 3, 4, 5, 6, 7]);
}

/** Generic mixed element type. Usefull only for its name. **/
proc Mixed_t()
{
    return (MIXED, "MIXED", -1, -1, -1, [ElementTypeNull], [-1], [-1], [ElementTypeNull], [-1], [-1]);
}

/** Returns the tuple definition of an element. Defined as Element_t(type:c_int, name:string, numNodes:int, numEdges:int,numFacets:int,edgeType:[0..#numEdges]c_int,ptrToEdge:[0..numEdges]int,edgeList:[0..#ptrToEdge[numEdges]]int,facetType:[0..#numFacets]c_int,ptrToFacet:[0..numFacets]int,facetList:[0..#ptrToFacet[numFacets]]int).

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Tuple definition of an element depending on its type: :type:`Node_t <ElementTopology.Node_t>`, :type:`Bar2_t <ElementTopology.Bar2_t>`, :type:`Tri3_t <ElementTopology.Tri3_t>`, :type:`Quad4_t <ElementTopology.Quad4_t>`, :type:`Tetra4_t <ElementTopology.Tetra4_t>`, :type:`Pyra5_t <ElementTopology.Pyra5_t>`, :type:`Penta6_t <ElementTopology.Penta6_t>`, :type:`Hexa8_t <ElementTopology.Hexa8_t>`.
    :rtype: (c_int, string, int, int, int, [c_int], [int], [int], [c_int], [int], [int])
 **/
proc getElementDefinition(elemType : c_int)
{
    select elemType
    {
        when NODE do return Node_t();
        when BAR_2 do return Bar2_t();
        when TRI_3 do return Tri3_t();
        when QUAD_4 do return Quad4_t();
        when TETRA_4 do return Tetra4_t();
        when PYRA_5 do return Pyra5_t();
        when PENTA_6 do return Penta6_t();
        when HEXA_8 do return Hexa8_t();
        when MIXED do return Mixed_t();
        otherwise do halt("Unknown element : ", elemType);
    }

    // Should never return here
    return (ElementTypeNull, "ElementTypeNull", -1, -1, -1, [ElementTypeNull], [-1], [-1], [ElementTypeNull], [-1], [-1]);
}

/** Returns the name of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Name of a topological element.
    :rtype: string
**/
proc getElementName(elemType : c_int) : string
{
    return getElementDefinition(elemType)(1);
}

/** Returns the number of nodes of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Number of nodes of a topological element.
    :rtype: int
**/
proc getElementNumNodes(elemType : c_int) : int
{
    return getElementDefinition(elemType)(2);
}

/** Returns the number of edges of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Number of edges of a topological element.
    :rtype: int
**/
proc getElementNumEdges(elemType : c_int) : int
{
    return getElementDefinition(elemType)(3);
}

/** Returns the number of facets of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Number of facets of a topological element.
    :rtype: int
**/
proc getElementNumFacets(elemType : c_int, dim : int = 3) : int
{
    if (dim == 2)
    {
        return getElementNumEdges(elemType);
    }
    else if (dim == 1)
    {
        return getElementNumNodes(elemType);
    }

    return getElementDefinition(elemType)(4);
}

/** Returns the list of edge types of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Array containing the edge types contained in the element. Size numEdges. Domain [0..#numEdges]
    :rtype: [0..#numEdges] int
**/
proc getElementEdgeType(elemType : c_int)
{
    return getElementDefinition(elemType)(5);
}

/** Returns the list of starting index for each edge in the edge-to-node connectivity of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Array containing the starting index for each edge in the edge-to-node connectivity with the last element containing the length of the edge-to-node connectivity array. Size numEdges+1. Domain [0..numEdges]
    :rtype: [0..numEdges] int
**/
proc getElementPtrToEdge(elemType : c_int)
{
    return getElementDefinition(elemType)(6);
}

/** Returns the list of edge-to-node connectivity of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Array containing the list of edge-to-node connectivity. Size ptrToEdge[numEdges]. Domain [0..#ptrToEdge[numEdges]]
    :rtype: [0..#ptrToEdge[numEdges]] int
**/
proc getElementEdgeList(elemType : c_int)
{
    return getElementDefinition(elemType)(7);
}

/** Returns the list of facet types of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Array containing the facet types contained in the element. Size numFacets. Domain [0..#numFacets]
    :rtype: [0..#numFacets] int
**/
proc getElementFacetType(elemType : c_int, dim : int = 3)
{
    if (dim == 2)
    {
        return getElementEdgeType(elemType);
    }
    else if (dim == 1)
    {
        if (elemType == BAR_2)
        {
            return [NODE, NODE];
        }
        else if (elemType == NODE)
        {
            return [NODE];
        }
        else
        {
            halt("Wrong element type of dim 1: ", elemType);
            return [ElementTypeNull];
        }
    }

    return getElementDefinition(elemType)(8);
}

/** Returns the list of starting index for each facet in the facet-to-node connectivity of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Array containing the starting index for each facet in the facet-to-node connectivity with the last element containing the length of the facet-to-node connectivity array. Size numFacets+1. Domain [0..numFacets]
    :rtype: [0..numFacets] int
**/
proc getElementPtrToFacet(elemType : c_int, dim : int = 3)
{
    if (dim == 2)
    {
        return getElementPtrToEdge(elemType);
    }
    else if (dim == 1)
    {
        if (elemType == BAR_2)
        {
            return [0, 1, 2];
        }
        else if (elemType == NODE)
        {
            return [0, 1];
        }
        else
        {
            halt("Wrong element type of dim 1: ", elemType);
            return [-1];
        }
    }

    return getElementDefinition(elemType)(9);
}

/** Returns the list of facet-to-node connectivity of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Array containing the list of facet-to-node connectivity. Size ptrToFacet[numFacets]. Domain [0..#ptrToFacet[numFacets]]
    :rtype: [0..#ptrToFacet[numFacets]] int
**/
proc getElementFacetList(elemType : c_int, dim : int = 3)
{
    if (dim == 2)
    {
        return getElementEdgeList(elemType);
    }
    else if (dim == 1)
    {
        if (elemType == BAR_2)
        {
            return [0, 1];
        }
        else if (elemType == NODE)
        {
            return [0];
        }
        else
        {
            halt("Wrong element type of dim 1: ", elemType);
            return [-1];
        }
    }

    return getElementDefinition(elemType)(10);
}

/** Find the topological dimension of a topological element

    :arg elemType: Integer code for element type from CGNS standard.
    :type elemType: c_int

    :returns: Topological dimension of a topological element.
    :rtype: int
**/
proc getElementDimension(elemType : c_int) : int
{
    var elem = getElementDefinition(elemType);

    if (elem(4) > 1) then return 3;
    else if (elem(4) == 1) then return 2;
    else if (elem(3) == 1) then return 1;
    else return 0;
}
} // module ElementTopology
