package tiler.tiling;

import java.util.Arrays;
import java.util.List;

import javafx.geometry.Point3D;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Rotate;
import tiler.core.dsymbols.Geometry;

public class Circle3D {

	public static Point3D[] circle(Point3D center, Point3D orientation, double radius, int N, Geometry geom) {

		Point3D[] coordinates = new Point3D[N];

		if (geom == Geometry.Hyperbolic) {
			coordinates = Tools.hyperbolicCircleCoordinates(center, orientation, radius, N);
		} else {

			// finds coordinates of a regular n-sided polygon in the x y plane with center
			// at 0
			for (int n = 0; n < N; n++) {
				coordinates[n] = new Point3D(radius * Math.cos(2 * Math.PI * n / N),
						radius * Math.sin(2 * Math.PI * n / N), 0);
			}

			// Finds normal vector
			Point3D normal = Tools.getNormalVector(center, geom); // new z Axis
			orientation = orientation.normalize(); // new x Axis
			Point3D newYAxis = normal.crossProduct(orientation).normalize(); // new y Axis
			// Transform points with Matrix multiplication. Affine Transformation to the
			// plane of the center point
			for (int n = 0; n < N; n++) {
				double newX = orientation.getX() * coordinates[n].getX() + newYAxis.getX() * coordinates[n].getY()
						+ normal.getX() * coordinates[n].getZ();
				double newY = orientation.getY() * coordinates[n].getX() + newYAxis.getY() * coordinates[n].getY()
						+ normal.getY() * coordinates[n].getZ();
				double newZ = orientation.getZ() * coordinates[n].getX() + newYAxis.getZ() * coordinates[n].getY()
						+ normal.getZ() * coordinates[n].getZ();
				coordinates[n] = new Point3D(newX, newY, newZ).add(center);
			}
		}
		return coordinates;

	}

	public static TriangleMesh CircleMesh(Point3D center, Point3D[] coordinates, Geometry geom) {

		Point3D[] points3d = new Point3D[coordinates.length + 1];
		points3d[0] = center;
		for (int i = 0; i < coordinates.length; i++) {
			points3d[i + 1] = coordinates[i];
		}

		int N = coordinates.length;
		int[] fac = new int[6 * N];
		int counter = 2;
		for (int i = 0; i < 6 * N; i = i + 6) {
			if (counter != N + 1) {
				fac[i] = 0;
				fac[i + 1] = 0;
				fac[i + 2] = counter - 1;
				fac[i + 3] = 1;
				fac[i + 4] = counter;
				fac[i + 5] = 2;
				counter++;
			} else {
				fac[i] = 0;
				fac[i + 1] = 0;
				fac[i + 2] = counter - 1;
				fac[i + 3] = 1;
				fac[i + 4] = 1;
				fac[i + 5] = 2;

			}
		}

		if (geom != Geometry.Spherical) {
			fac = tiler.tiling.FundamentalDomain.invertOrientation(fac);
		}

		float[] points = new float[3 * points3d.length];

		for (int i = 0; i < points3d.length; i++) {
			points[3 * i] = (float) points3d[i].getX();
			points[3 * i + 1] = (float) points3d[i].getY();
			points[3 * i + 2] = (float) points3d[i].getZ();
		}

		final float[] texCoords = { 0.5f, 0, 0, 0, 1, 1 };
		int[] smoothing = new int[fac.length / 6];
		Arrays.fill(smoothing, 1);

		TriangleMesh mesh = new TriangleMesh();
		mesh.getPoints().addAll(points);
		mesh.getTexCoords().addAll(texCoords);
		mesh.getFaces().addAll(fac);
		mesh.getFaceSmoothingGroups().addAll(smoothing);

		return mesh;

	}

	public static TriangleMesh CircleMesh(Point3D center, Point3D[] coordinates, Geometry geom, double above) {

		Point3D[] points3d = new Point3D[coordinates.length + 1];
		points3d[0] = center;
		for (int i = 0; i < coordinates.length; i++) {
			points3d[i + 1] = coordinates[i];
		}

		int N = coordinates.length;
		int[] fac = new int[6 * N];
		int counter = 2;
		for (int i = 0; i < 6 * N; i = i + 6) {
			if (counter != N + 1) {
				fac[i] = 0;
				fac[i + 1] = 0;
				fac[i + 2] = counter - 1;
				fac[i + 3] = 1;
				fac[i + 4] = counter;
				fac[i + 5] = 2;
				counter++;
			} else {
				fac[i] = 0;
				fac[i + 1] = 0;
				fac[i + 2] = counter - 1;
				fac[i + 3] = 1;
				fac[i + 4] = 1;
				fac[i + 5] = 2;

			}
		}

		for (int i = 0; i < points3d.length; i++) {
			Point3D normal = tiler.tiling.Tools.getNormalVector(points3d[i], geom);
			points3d[i] = points3d[i].add(normal.multiply(above));
		}

		if (geom != Geometry.Spherical) {
			fac = tiler.tiling.FundamentalDomain.invertOrientation(fac);
		}

		float[] points = new float[3 * points3d.length];

		for (int i = 0; i < points3d.length; i++) {
			points[3 * i] = (float) points3d[i].getX();
			points[3 * i + 1] = (float) points3d[i].getY();
			points[3 * i + 2] = (float) points3d[i].getZ();
		}

		final float[] texCoords = { 0.5f, 0, 0, 0, 1, 1 };
		int[] smoothing = new int[fac.length / 6];
		Arrays.fill(smoothing, 1);

		TriangleMesh mesh = new TriangleMesh();
		mesh.getPoints().addAll(points);
		mesh.getTexCoords().addAll(texCoords);
		mesh.getFaces().addAll(fac);
		mesh.getFaceSmoothingGroups().addAll(smoothing);

		return mesh;

	}

	public static MeshView CircleMeshView(TriangleMesh mesh, Color color) {
		MeshView meshView = new MeshView(mesh);
		PhongMaterial material = new PhongMaterial(color);
		meshView.setMaterial(material);
		return meshView;
	}

}
