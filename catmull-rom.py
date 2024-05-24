"""
catmull-rom.py
"""

import numpy as np

__all__ = ["CatmullRom"]


class CatmullRom :

################################################################################
    def num_segments(points: tuple) -> int:
        """
        Returns the amount of segments included in the point chain. There are 
        n-1 segments for n points, and 2 of them are "ghost segments", so we
        subtract 3.
        :param points: List of control points.
        :return: The number of segments in the spline.
        """
        
        return len(points) - 3

    def tj(ti: float, pi: tuple, pj: tuple, alpha:float) -> float:
        """
        Returns the value of the next t depending on the distance between current
        p and next p, as well as some alpha to determine knot parameterization.
        :param ti: The current t.
        :param pi: The current control point.
        :param pj: The next control point.
        :param alpha: 0.0 for uniform spline, 0.5 for centripetal spline, 1.0
            for chordal spline.
        :return: The value of the next t.
        """
        # felt cute might rewrite
        return ((((pi[0] - pj[0]) ** 2) + ((pi[1] - pj[1]) ** 2) + ((pi[2] - pj[2]) ** 2)) ** 0.5) ** alpha + ti
    
    def spline(
        p0: tuple, 
        p1: tuple, 
        p2: tuple, 
        p3: tuple, 
        alpha: float, 
        num_points: int
    ) -> list[tuple[float, float, float]]:
        """
        Returns points that make up the spline segment.
        :param p0, p1, p2, p3: The points that define the Catmull-Rom spline.
        :param alpha: 0.0 for uniform spline, 0.5 for centripetal spline, 1.0
            for chordal spline.
        :param num_points: Number of points that should make up the spline.
        :return: List of points that make up the spline
        """
        t0 : float = 0.0
        t1 : float = tj(t0, p0, p1)
        t2 : float = tj(t1, p1, p2)
        t3 : float = tj(t2, p2, p3)

        
