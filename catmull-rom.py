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

    def tj(ti: float, Pi: tuple, Pj: tuple, alpha:float) -> float:
        """
        Returns the value of the next t depending on the distance between current
        p and next p, as well as some alpha to determine knot parameterization.
        :param ti: The current t.
        :param Pi: The current control point.
        :param Pj: The next control point.
        :param alpha: 0.0 for uniform spline, 0.5 for centripetal spline, 1.0
            for chordal spline.
        :return: The value of the next t.
        """
        # felt cute might rewrite
        return ((((Pi[0] - Pj[0]) ** 2) + ((Pi[1] - Pj[1]) ** 2) + ((Pi[2] - Pj[2]) ** 2)) ** 0.5) ** alpha + ti
    
    def spline(
        self,
        P0: tuple, 
        P1: tuple, 
        P2: tuple, 
        P3: tuple, 
        alpha: float, 
        num_points: int
    ) -> list[tuple[float, float, float]]:
        """
        Returns points that make up the spline segment.
        :param P0, P1, P2, P3: The points that define the Catmull-Rom spline.
        :param alpha: 0.0 for uniform spline, 0.5 for centripetal spline, 1.0
            for chordal spline.
        :param num_points: Number of points that should make up the spline.
        :return: List of points that make up the spline
        """
        t0 : float = 0.0
        t1 : float = self.tj(t0, P0, P1)
        t2 : float = self.tj(t1, P1, P2)
        t3 : float = self.tj(t2, P2, P3)
        t = np.linspace(t1, t2, num_points).reshape(num_points, 1)

        A1 = ((t1 - t) / (t1 - t0)) * P0 + ((t - t0) / (t1 - t0)) * P1
        A2 = ((t2 - t) / (t2 - t1)) * P1 + ((t - t1) / (t2 - t0)) * P2
        A3 = ((t3 - t) / (t3 - t2)) * P2 + ((t - t2) / (t3 - t0)) * P3

        B1 = ((t2 - t) / (t2 - t0)) * A1 + ((t - t0) / (t2 - t0)) * A2
        B2 = ((t3 - t) / (t3 - t1)) * A2 + ((t - t1) / (t3 - t1)) * A3

        C = ((t2 - t) / (t2 - t1)) * B1 + ((t - t1) / (t2 - t1)) * B2

        return C