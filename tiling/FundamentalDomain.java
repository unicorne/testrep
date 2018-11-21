package tiler.tiling;

import javafx.geometry.Point2D;
import javafx.geometry.Point3D;
import javafx.scene.AmbientLight;
import javafx.scene.Group;
import javafx.scene.PerspectiveCamera;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.Circle;
import javafx.scene.shape.CullFace;
import javafx.scene.shape.Line;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.Sphere;
import javafx.scene.shape.Shape;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.scene.transform.Affine;
import javafx.scene.transform.Transform;
import javafx.scene.transform.Translate;
import tiler.core.dsymbols.DSymbol;
import tiler.core.dsymbols.FDomain;
import tiler.core.dsymbols.Geometry;
import tiler.core.fundamental.utils.WrapInt;
import tiler.core.fundamental.utils.Gaps;

import tiler.core.fundamental.tools.Line3D;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.Random;
import java.util.List;
import java.util.ArrayList;

import com.sun.javafx.css.CalculatedValue;

/**
 * builds fundamental domain in JavaFX Created by huson on 4/5/16.
 */
public class FundamentalDomain {
	/**
	 * construct a fundamental domain
	 *
	 * @param dsymbol
	 *            Delaney symbol from old DH code
	 * @param fDomain
	 *            domain computed by KW
	 * @return fundamental domain
	 */

	public static Group buildFundamentalDomain(final DSymbol dsymbol, final FDomain fDomain) {

		final Group group = new Group();

		final Color[] colors = new Color[fDomain.size() + 1];

		final BitSet set = new BitSet();
		final BitSet visit = new BitSet(fDomain.size());
		final Random random = new Random(666);
		// set colors
		for (int a = 1; a <= dsymbol.size(); a = dsymbol.nextOrbit(0, 1, a, set)) {
			final Color color = new Color(random.nextDouble(), random.nextDouble(), random.nextDouble(), 1);
			// final Color color = new Color(1, 1, 1, 1); //all white
			dsymbol.visitOrbit(0, 1, a, new DSymbol.OrbitVisitor() {
				public void visit(int a) {
					colors[a] = color;
				}
			});
		}

		// For lines
		float[][] pointsmatrix; // saves points for mesh drawing for each fDomain
		Point3D[][] points3dmatrix; // saves point3ds for mesh drawing for each fDomain
		int[][] linepointsmatrix;// saves the position of points that are on a line that needs to be drawn
		double linesize = 2;
		double edgesize = 2;
		int edgefine = 8; // defines how smooth the edges are
		final Color linecolor = Color.WHITE;
		final Color lineedgecolor = Color.BLACK;
		final Color edgecolor = Color.WHITE;
		double linesAbove = 0.5; // defines the height of the line above the faces
		// Booleans
		boolean drawFaces = true;
		boolean drawLines = true;
		boolean drawEdges = true;
		boolean drawPoints = false; // draws points on 2 line

		boolean drawLinesWithEdges = false; // Attempt to draw lines with edges

		boolean drawEucledianEdges = true; // draws edges with circles
		boolean drawEucledianEdgesSphere = false; // draws edges with sphere objects

		boolean drawSphericalEdges = true;
		boolean drawSphericalEdgesFine = false; // more edges less performance
		boolean drawSphericalEdgeSphere = false;
		boolean drawSphericalEdgesFineSphere = false;

		boolean drawHyperbolicEdges = false;
		boolean drawHyperbolicEdgesFine = true;
		boolean drawHyperbolicEdgesSphere = false;
		boolean drawHyperbolicEdgesSphereFine = false;

		boolean HyperbolicAccurateLineSize = false; // hyperbolic linesize is adjusted

		// defines length of the matrix
		Geometry geom = fDomain.getGeometry();
		if (geom == Geometry.Spherical) {
			points3dmatrix = new Point3D[fDomain.size()][1026]; // number of points that create each spherical fDomain
			linepointsmatrix = new int[fDomain.size()][33]; // number of points on the line that needs to be drawn
		} else if (geom == Geometry.Euclidean) {
			points3dmatrix = new Point3D[fDomain.size()][7];
			linepointsmatrix = new int[fDomain.size()][3];
		} else {
			points3dmatrix = new Point3D[fDomain.size()][13];
			linepointsmatrix = new int[fDomain.size()][9];
		}

		// construct triangles as meshes:

		final int orientation = (computeWindingNumber(fDomain.getVertex3D(0, 1), fDomain.getVertex3D(1, 1),
				fDomain.getVertex3D(2, 1)) < 0 ? fDomain.getOrientation(1) : -fDomain.getOrientation(1));

		for (int a = 1; a <= fDomain.size(); a++) {
			final float[] points;
			final Point3D[] points3d; // points that create the triangles

			final int[] faces;
			final int[] fac;

			final int[] smoothing;

			final int[] linepoints;

			if (geom == Geometry.Spherical) {

				// Spherical

				int depth = 4; // 4^5 = 1024

				fac = new int[(int) Math.pow(4, (depth + 1)) * 6];
				points3d = new Point3D[1026]; // 3, 6, 66, 258, 1026 // size of points array dependent on depth

				WrapInt p = new WrapInt(0);
				WrapInt f = new WrapInt(0);

				points3d[p.incrementInt()] = fDomain.getVertex3D(0, a);
				points3d[p.incrementInt()] = fDomain.getVertex3D(1, a);
				points3d[p.incrementInt()] = fDomain.getVertex3D(2, a);
				points3d[p.incrementInt()] = fDomain.getEdgeCenter3D(0, a);
				points3d[p.incrementInt()] = fDomain.getEdgeCenter3D(1, a);
				points3d[p.incrementInt()] = fDomain.getEdgeCenter3D(2, a);

				// Iterative triangle mesh generator
				class triangle {

					private boolean orientationUp;
					private int pointA, pointB, pointC;
					private int depth;
					private triangle tri1;
					private triangle tri2;
					private triangle tri3;
					private triangle tri4;

					triangle(boolean orientationUp, int pointA, int pointB, int pointC, int depth) {
						this.orientationUp = orientationUp;
						this.pointA = pointA;
						this.pointB = pointB;
						this.pointC = pointC;
						this.depth = depth;

						if (this.depth > 0) {
							int midAB = p.incrementInt();
							points3d[midAB] = Tools.midpoint3D(geom, points3d[pointA], points3d[pointB]);
							int midAC = p.incrementInt();
							points3d[midAC] = Tools.midpoint3D(geom, points3d[pointA], points3d[pointC]);
							int midBC = p.incrementInt();
							points3d[midBC] = Tools.midpoint3D(geom, points3d[pointB], points3d[pointC]);

							this.tri1 = new triangle(this.orientationUp, this.pointA, midAB, midAC, --this.depth);
							this.tri2 = new triangle(this.orientationUp, midAB, this.pointB, midBC, this.depth);
							this.tri3 = new triangle(this.orientationUp, midAC, midBC, this.pointC, this.depth);

							if (this.orientationUp) {
								this.tri4 = new triangle(!this.orientationUp, midAB, midBC, midAC, this.depth);
							} else {
								this.tri4 = new triangle(!this.orientationUp, midAC, midAB, midBC, this.depth);
							}
						} else {
							int facPos = 6 * f.incrementInt();
							fac[facPos] = pointA;
							fac[facPos + 1] = 0;
							fac[facPos + 2] = pointB;
							fac[facPos + 3] = 1;
							fac[facPos + 4] = pointC;
							fac[facPos + 5] = 2;
						}
					}
				}

				// clockwise orientation
				new triangle(true, 0, 4, 5, depth);
				new triangle(true, 5, 3, 1, depth);
				new triangle(true, 4, 2, 3, depth);
				new triangle(false, 4, 3, 5, depth);

				int[] pointsOf2Edge = { 0, 1, 5, 7, 10, 13, 16, 22, 43, 46, 52, 136, 139, 142, 148, 169, 172, 178, 262,
						265, 268, 271, 277, 298, 301, 307, 391, 394, 397, 403, 424, 427, 433 }; // getPositionArray(bolArray);

				int[] pointsOf2EdgeSorted = { 0, 16, 13, 22, 10, 46, 43, 52, 7, 142, 139, 148, 136, 172, 169, 178, 5,
						271, 268, 277, 265, 301, 298, 307, 262, 397, 394, 403, 391, 427, 424, 433, 1 };

				linepoints = pointsOf2EdgeSorted;
				// double[][] linepointsdistance = new double[1026][2];
				// double[][] linepointsdistance2 = new double[1026][2];
				//
				// double[][] linepointsdistance3 ={{0.0, 0.0},
				// {1.0, 9.219244256272484},
				// {5.0, 4.609622128136366},
				// {7.0, 2.3048110640685406},
				// {10.0, 1.1524055320345405},
				// {13.0, 0.5762027660157822},
				// {16.0, 0.288101383005441},
				// {22.0, 0.8643041490261935},
				// {43.0, 1.72860829805166},
				// {46.0, 1.4405069150418683},
				// {52.0, 2.016709681059741},
				// {136.0, 3.457216596102366},
				// {139.0, 2.8810138300846284},
				// {142.0, 2.592912447077225},
				// {148.0, 3.169115213094037},
				// {169.0, 4.033419362120038},
				// {172.0, 3.7453179791108875},
				// {178.0, 4.321520745127799},
				// {262.0, 6.9144331922046876},
				// {265.0, 5.762027660170338},
				// {268.0, 5.185824894153019},
				// {271.0, 4.897723511145091},
				// {277.0, 5.473926277162035},
				// {298.0, 6.338230426187096},
				// {301.0, 6.050129043179344},
				// {307.0, 6.626331809196125},
				// {391.0, 8.06683872423854},
				// {394.0, 7.4906359582212065},
				// {397.0, 7.2025345752132095},
				// {403.0, 7.778737341230291},
				// {424.0, 8.643041490255458},
				// {427.0, 8.354940107246374},
				// {433.0, 8.931142873264196}};
				//
				//
				// for(int i=0;i<33;i++) {
				// linepointsdistance[i][0] = positions2[i];
				// linepointsdistance[i][1] =
				// sphericalDistance(points3d[0],points3d[positions2[i]],100);
				// }
				//
				// for(int i = 0;i<33;i++) {
				// System.out.println("" + Arrays.toString(linepointsdistance3[i]));
				// }
				//
				// linepointsdistance3 = sort2DArrayBasedOnSecondColumn(linepointsdistance3);
				//
				// for(int i = 0;i<33;i++) {
				// System.out.println("sorted" + Arrays.toString(linepointsdistance3[i]));
				// }
				//

				// outsource
				linepointsmatrix[a - 1] = linepoints;
				points3dmatrix[a - 1] = points3d;

				// tries to get points
				if (false) {

					int[] BooleanPoints = new int[1026];
					Point3D start = points3d[0];
					Point3D end = points3d[1];
					double startenddist = start.distance(end);
					double dist11 = start.distance(points3d[5]);
					double dist22 = end.distance(points3d[5]);
					double dist33 = dist11 + dist22;
					double diff = Math.abs(dist33 - startenddist);
					BooleanPoints[1] = 1;
					BooleanPoints[0] = 1;

					for (int i = 2; i < 1026; i++) {
						double dist1 = start.distance(points3d[i]);
						double dist2 = end.distance(points3d[i]);
						double dist3 = dist1 + dist2;

						if (Math.abs(startenddist - dist3) <= diff) {
							BooleanPoints[i] = 1;
						} else {
							BooleanPoints[i] = 0;
						}
					}
					int counter = 0;
					for (int i = 0; i < 1026; i++) {
						if (BooleanPoints[i] == 1) {
							counter++;
						}
					}

					// System.out.println("counter: " + counter);
					int[] positions = new int[counter];
					int j = 0;
					for (int i = 0; i < 1026; i++) {
						if (BooleanPoints[i] == 1) {
							positions[j] = i;
							j++;
						}
					}
					// System.out.println(Arrays.toString(positions));

				}

			} else if (geom == Geometry.Euclidean) {

				// Euclidean

				/// Original mesh structure
				points3d = new Point3D[7];
				int p = 0;
				for (int i = 0; i <= 2; i++) {
					points3d[p++] = fDomain.getVertex3D(i, a);
				}
				for (int i = 0; i <= 2; i++) {
					points3d[p++] = fDomain.getEdgeCenter3D(i, a);
				}
				points3d[p++] = fDomain.getChamberCenter3D(a);
				int[] original = new int[] { 0, 0, 6, 1, 5, 2, // v0 cc e2
						1, 0, 5, 1, 6, 2, // v1 e2 cc
						1, 0, 6, 1, 3, 2, // v1 cc e0
						2, 0, 3, 0, 6, 2, // v2 e0 cc
						2, 0, 6, 1, 4, 2, // v2 cc e1
						0, 0, 4, 1, 6, 2 // v0 e1 cc
				};

				// Reduced mesh structure: Mesh consists only of 2 triangles
				// points3d = new Point3D[4]; //4
				//
				// int p = 0;
				// for (int i = 0; i <= 2; i++) {
				// points3d[p++] = fDomain.getVertex3D(i,a);
				// }
				// points3d[p++] = fDomain.getEdgeCenter3D(2,a);
				//
				// int[] original = new int[]{
				// 0, 0, 2, 1, 3, 2, //v0 v2 e2
				// 2, 0, 1, 1, 3, 2, //v2 v1 e2
				// };

				fac = original;

				// outsource
				int[] positions = { 0, 5, 1 };
				linepoints = positions;
				points3dmatrix[a - 1] = points3d;
				linepointsmatrix[a - 1] = linepoints;

			} else {

				// hyperbolic

				points3d = new Point3D[13];

				int p = 0;

				for (int i = 0; i <= 2; i++) {
					points3d[p++] = fDomain.getVertex3D(i, a);
				}
				for (int i = 0; i <= 2; i++) {
					points3d[p++] = fDomain.getEdgeCenter3D(i, a);
				}
				points3d[p++] = fDomain.getChamberCenter3D(a);

				// hyper
				points3d[p++] = Tools.midpoint3D(geom, points3d[0], points3d[5]);
				points3d[p++] = Tools.midpoint3D(geom, points3d[5], points3d[1]);
				points3d[p++] = Tools.midpoint3D(geom, points3d[0], points3d[7]);
				points3d[p++] = Tools.midpoint3D(geom, points3d[7], points3d[5]);
				points3d[p++] = Tools.midpoint3D(geom, points3d[5], points3d[8]);
				points3d[p++] = Tools.midpoint3D(geom, points3d[8], points3d[1]);

				int[] hyper = new int[] { 0, 0, 6, 1, 9, 2, //
						9, 0, 6, 1, 7, 2, //
						7, 0, 6, 1, 10, 2, //
						10, 0, 6, 1, 5, 2, //
						5, 0, 6, 1, 11, 2, //
						11, 0, 6, 1, 8, 2, //
						8, 0, 6, 1, 12, 2, //
						12, 0, 6, 1, 1, 2, //
						0, 0, 4, 1, 6, 2, //
						4, 0, 2, 1, 6, 2, //
						2, 0, 3, 1, 6, 2, //
						6, 0, 3, 1, 1, 2 //
				};

				fac = hyper;

				// outsource
				int[] pointsOf2EdgeSorted = { 0, 9, 7, 10, 5, 11, 8, 12, 1 };
				linepoints = pointsOf2EdgeSorted;
				points3dmatrix[a - 1] = points3d;
				linepointsmatrix[a - 1] = linepoints;

			} // end of geometric cases

			points = new float[3 * points3d.length];

			for (int i = 0; i < points3d.length; i++) {
				points[3 * i] = (float) points3d[i].getX();
				points[3 * i + 1] = (float) points3d[i].getY();
				points[3 * i + 2] = (float) points3d[i].getZ();
			}

			if (fDomain.getOrientation(a) == orientation) {
				faces = fac;
			} else {
				faces = invertOrientation(fac);
			}

			smoothing = new int[faces.length / 6];
			Arrays.fill(smoothing, 1);

			final float[] texCoords = { 0.5f, 0, 0, 0, 1, 1 };

			// outsource
			pointsmatrix = new float[fDomain.size()][3 * points3d.length];
			pointsmatrix[a - 1] = points;

			// Draw Faces
			if (drawFaces) {
				TriangleMesh mesh = new TriangleMesh();
				mesh.getPoints().addAll(points);
				mesh.getTexCoords().addAll(texCoords);
				mesh.getFaces().addAll(faces);
				mesh.getFaceSmoothingGroups().addAll(smoothing);
				MeshView meshView = new MeshView(mesh);
				meshView.setMesh(mesh);
				PhongMaterial material = new PhongMaterial(colors[a]);
				// material.setSpecularColor(Color.YELLOW);
				meshView.setMaterial(material);
				group.getChildren().addAll(meshView);
			}

			// Draw Points
			if (drawPoints) {
				for (int i = 0; i < linepoints.length; i++) {
					Sphere pdot = new Sphere(1);
					pdot.setMaterial(new PhongMaterial(Color.BLACK));
					pdot.getTransforms().add(new Translate(points3d[linepoints[i]].getX(),
							points3d[linepoints[i]].getY(), points3d[linepoints[i]].getZ()));
					group.getChildren().add(pdot);
				}
			}
		}

		// Draw Lines
		if (true) {
			for (int a = 1; a <= fDomain.size(); a++) {
				// draw lines out
				if (drawLines) {
					PhongMaterial linematerial = new PhongMaterial(linecolor);
					PhongMaterial edgematerial = new PhongMaterial(lineedgecolor);

					if (HyperbolicAccurateLineSize && fDomain.getGeometry() == Geometry.Hyperbolic) {
						Point3D refPoint = fDomain.getChamberCenter3D(1).multiply(0.01);
						Point3D origin = new Point3D(0, 0, 1);
						double w = 0.01;
						double h = (1 + w * w) / (1 - w * w);
						// Length of translation
						double t = Tools.distance(fDomain, refPoint, origin);
						// Affine translation:
						Affine translateT = new Affine(Math.cosh(t), Math.sinh(t), 0, Math.sinh(t), Math.cosh(t), 0); // Translation
																														// along
																														// x-axis
						Point2D x = translateT.transform(0, 1);
						Point2D y = translateT.transform((1 + h) * w, h);

						linesize = 100 * (y.getX() / (1 + y.getY()) - x.getX() / (1 + x.getY()));
						edgesize = linesize;
					}

					if (!visit.get(a)) {
						visit.set(a);
						//visit.set(dsymbol.getS2(a));
						for (int i = 0; i < linepointsmatrix[a - 1].length - 1; i++) {
							int k = linepointsmatrix[a - 1][i];
							int j = linepointsmatrix[a - 1][i + 1];
							if (drawLinesWithEdges) {
								edgesize = linesize * 0.6;
								TriangleMesh linemesh = new TriangleMesh();
								TriangleMesh linemeshedges = new TriangleMesh();
								linemesh = tiler.core.fundamental.tools.Line3D.connectEdges(points3dmatrix[a - 1][k],
										points3dmatrix[a - 1][j], fDomain, linesize, 0.2, linesAbove)[0];
								linemeshedges = tiler.core.fundamental.tools.Line3D.connectEdges(
										points3dmatrix[a - 1][k], points3dmatrix[a - 1][j], fDomain, linesize, 0.2,
										linesAbove)[1];

								MeshView lineView = new MeshView();
								lineView.setMesh(linemesh);
								lineView.setMaterial(linematerial);

								MeshView lineView2 = new MeshView();
								lineView2.setMesh(linemeshedges);
								lineView2.setMaterial(edgematerial);

								group.getChildren().addAll(lineView);
								group.getChildren().addAll(lineView2);
							} else {
								TriangleMesh linemesh = new TriangleMesh();
								linemesh = tiler.core.fundamental.tools.Line3D.connect(points3dmatrix[a - 1][k],
										points3dmatrix[a - 1][j], fDomain, linesize, linesAbove);

								MeshView lineView = new MeshView();
								lineView.setMesh(linemesh);
								lineView.setMaterial(linematerial);

								group.getChildren().addAll(lineView);
							}

						}
					}
				}
				// handle edges
				if (drawEdges) {
					if (geom == Geometry.Euclidean) {
						if (drawEucledianEdges) {
							int[] EdgesToDraw = { 0, 5, 1 };
							for (int i = 0; i < EdgesToDraw.length; i++) {
								Point3D center = points3dmatrix[a - 1][EdgesToDraw[i]].add(0, 0, 1);
								Point3D diff = new Point3D(1, 0, 0); // points3dmatrix[a-1][5].add(0,0,1).subtract(center);
								Point3D[] coordinates = tiler.core.fundamental.tools.Circle3D.circle(center, diff,
										edgesize, edgefine, geom);
								TriangleMesh mesh = tiler.core.fundamental.tools.Circle3D.CircleMesh(center,
										coordinates, geom);
								MeshView meshView = tiler.core.fundamental.tools.Circle3D.CircleMeshView(mesh,
										edgecolor);
								group.getChildren().addAll(meshView);
							}

						} else if (drawEucledianEdgesSphere) {
							int[] EdgesToDraw = {0, 5, 1};
							for (int i = 0; i < EdgesToDraw.length; i++) {
								Sphere pdot = new Sphere(1);
								pdot.setMaterial(new PhongMaterial(edgecolor));
								pdot.getTransforms().add(new Translate(points3dmatrix[a - 1][EdgesToDraw[i]].getX(),
										points3dmatrix[a - 1][EdgesToDraw[i]].getY(), 0));
								group.getChildren().add(pdot);
							}
						}
					}
					if (geom == Geometry.Spherical) {

						if (drawSphericalEdges) {
							int[] EdgesToDraw = { 0, 9, 16, 24,32};
							for (int i = 0; i < EdgesToDraw.length; i++) {
								int k;
								int l = EdgesToDraw[i];
								if (i == EdgesToDraw.length - 1) {
									k = EdgesToDraw[i - 1];
								} else {
									k = EdgesToDraw[i + 1];
								}
								Point3D normalVector = getSphericalNormal(
										points3dmatrix[a - 1][linepointsmatrix[a - 1][l]]);
								Point3D center = points3dmatrix[a - 1][linepointsmatrix[a - 1][l]].add(normalVector.multiply(linesAbove));
								Point3D diff = points3dmatrix[a - 1][linepointsmatrix[a - 1][k]].add(normalVector.multiply(linesAbove))
										.subtract(center);
								Point3D[] coordinates = tiler.core.fundamental.tools.Circle3D.circle(center, diff,
										edgesize, edgefine, geom);
								TriangleMesh mesh = tiler.core.fundamental.tools.Circle3D.CircleMesh(center,
										coordinates, geom);
								MeshView meshView = tiler.core.fundamental.tools.Circle3D.CircleMeshView(mesh,
										edgecolor);
								group.getChildren().addAll(meshView);
							}
						} else if (drawSphericalEdgesFine) {

							for (int i = 0; i < linepointsmatrix[a - 1].length; i++) {
								int k;
								if (i == linepointsmatrix[a - 1].length - 1) {
									k = i - 1;
								} else {
									k = i + 1;
								}
								Point3D normalVector = getSphericalNormal(
										points3dmatrix[a - 1][linepointsmatrix[a - 1][i]]);
								Point3D center = points3dmatrix[a - 1][linepointsmatrix[a - 1][i]].add(normalVector.multiply(linesAbove));
								Point3D diff = points3dmatrix[a - 1][linepointsmatrix[a - 1][k]].add(normalVector.multiply(linesAbove))
										.subtract(center);
								Point3D[] coordinates = tiler.core.fundamental.tools.Circle3D.circle(center, diff,
										edgesize, edgefine, geom);
								TriangleMesh mesh = tiler.core.fundamental.tools.Circle3D.CircleMesh(center,
										coordinates, geom);
								MeshView meshView = tiler.core.fundamental.tools.Circle3D.CircleMeshView(mesh,
										edgecolor);
								group.getChildren().addAll(meshView);
							}
						}else if (drawSphericalEdgeSphere) {
							int[] EdgesToDraw = { 0, 9, 16, 24,32};
							for (int i = 0; i < EdgesToDraw.length; i++) {
								int k=linepointsmatrix[a-1][EdgesToDraw[i]];
							Sphere pdot = new Sphere(1);
							pdot.setMaterial(new PhongMaterial(edgecolor));
							pdot.getTransforms().add(new Translate(points3dmatrix[a - 1][k].getX(),
									points3dmatrix[a - 1][k].getY(), points3dmatrix[a - 1][k].getZ()));
							group.getChildren().add(pdot);
							}

						
						} else if (drawSphericalEdgesFineSphere) {

							for (int i = 0; i < linepointsmatrix[a - 1].length; i++) {
								Sphere pdot = new Sphere(1);
								pdot.setMaterial(new PhongMaterial(edgecolor));
								pdot.getTransforms()
										.add(new Translate(points3dmatrix[a - 1][linepointsmatrix[a - 1][i]].getX(),
												points3dmatrix[a - 1][linepointsmatrix[a - 1][i]].getY(),
												points3dmatrix[a - 1][linepointsmatrix[a - 1][i]].getZ()));
								group.getChildren().add(pdot);
							}
						} 
					}
					if (geom == Geometry.Hyperbolic) {
						if (drawHyperbolicEdgesSphereFine && geom == Geometry.Hyperbolic) {
							for (int i = 0; i < linepointsmatrix[a - 1].length - 1; i++) {

								double above = 0.5;
								Point3D point0 = points3dmatrix[a - 1][linepointsmatrix[a - 1][i]];
								Point3D point1 = points3dmatrix[a - 1][linepointsmatrix[a - 1][i + 1]];
								Point3D diff = point1.subtract(point0).normalize();
								Point3D normalToSphere = new Point3D(2 * point0.getX(), 2 * point0.getY(),
										-2 * point0.getZ());
								Point3D normal = diff.crossProduct(normalToSphere.normalize());
								Point3D point2 = point0.add(normalToSphere.normalize().multiply(above));

								Sphere pdot = new Sphere(1);
								pdot.setMaterial(new PhongMaterial(edgecolor));
								pdot.getTransforms().add(new Translate(point2.getX(), point2.getY(), point2.getZ()));
								group.getChildren().add(pdot);
							}
						} else if (drawHyperbolicEdgesFine && geom == Geometry.Hyperbolic) {
							for (int i = 0; i < linepointsmatrix[a - 1].length; i++) {
								int k;
								if (i == linepointsmatrix[a - 1].length - 1) {
									k = i - 1;
								} else {
									k = i + 1;
								}
								Point3D normalVector = getHyperbolicNormal(
										points3dmatrix[a - 1][linepointsmatrix[a - 1][i]]);
								Point3D center = points3dmatrix[a - 1][linepointsmatrix[a - 1][i]].add(normalVector);
								Point3D diff = points3dmatrix[a - 1][linepointsmatrix[a - 1][k]].add(normalVector)
										.subtract(center);
								Point3D[] coordinates = tiler.core.fundamental.tools.Circle3D.circle(center, diff,
										edgesize, edgefine, geom);
								TriangleMesh mesh = tiler.core.fundamental.tools.Circle3D.CircleMesh(center,
										coordinates, geom);
								MeshView meshView = tiler.core.fundamental.tools.Circle3D.CircleMeshView(mesh,
										edgecolor);
								group.getChildren().addAll(meshView);
							}

						}
						if (drawHyperbolicEdges && geom == Geometry.Hyperbolic) {
							Point3D normalVector = getHyperbolicNormal(
									points3dmatrix[a - 1][linepointsmatrix[a - 1][0]]);
							Point3D normalVector2 = getHyperbolicNormal(
									points3dmatrix[a - 1][linepointsmatrix[a - 1][5]]);
							Point3D normalVector3 = getHyperbolicNormal(
									points3dmatrix[a - 1][linepointsmatrix[a - 1][1]]);

							Point3D center = points3dmatrix[a - 1][0].add(normalVector);
							Point3D diff = points3dmatrix[a - 1][5].add(normalVector).subtract(center);
							Point3D[] coordinates = tiler.core.fundamental.tools.Circle3D.circle(center, diff, edgesize,
									edgefine, geom);
							TriangleMesh mesh = tiler.core.fundamental.tools.Circle3D.CircleMesh(center, coordinates,
									geom);
							MeshView meshView = tiler.core.fundamental.tools.Circle3D.CircleMeshView(mesh, edgecolor);
							group.getChildren().addAll(meshView);
							////////////
							Point3D center2 = points3dmatrix[a - 1][5].add(normalVector3);
							Point3D diff2 = points3dmatrix[a - 1][1].add(normalVector3).subtract(center2);
							Point3D[] coordinates2 = tiler.core.fundamental.tools.Circle3D.circle(center2, diff2,
									edgesize, edgefine, geom);
							TriangleMesh mesh2 = tiler.core.fundamental.tools.Circle3D.CircleMesh(center2, coordinates2,
									geom);
							MeshView meshView2 = tiler.core.fundamental.tools.Circle3D.CircleMeshView(mesh2, edgecolor);
							group.getChildren().addAll(meshView2);
							////////////
							Point3D center3 = points3dmatrix[a - 1][1].add(normalVector2);
							Point3D diff3 = points3dmatrix[a - 1][5].add(normalVector2).subtract(center3);
							Point3D[] coordinates3 = tiler.core.fundamental.tools.Circle3D.circle(center3, diff3,
									edgesize, edgefine, geom);
							TriangleMesh mesh3 = tiler.core.fundamental.tools.Circle3D.CircleMesh(center3, coordinates3,
									geom);
							MeshView meshView3 = tiler.core.fundamental.tools.Circle3D.CircleMeshView(mesh3, edgecolor);
							group.getChildren().addAll(meshView3);

						} else if (drawHyperbolicEdgesSphere && geom == Geometry.Hyperbolic) {
							Sphere pdot = new Sphere(1);
							pdot.setMaterial(new PhongMaterial(edgecolor));
							pdot.getTransforms().add(new Translate(points3dmatrix[a - 1][0].getX(),
									points3dmatrix[a - 1][0].getY(), points3dmatrix[a - 1][0].getZ()));
							group.getChildren().add(pdot);

							Sphere pdot2 = new Sphere(1);
							pdot2.setMaterial(new PhongMaterial(edgecolor));
							pdot2.getTransforms().add(new Translate(points3dmatrix[a - 1][1].getX(),
									points3dmatrix[a - 1][1].getY(), points3dmatrix[a - 1][1].getZ()));
							group.getChildren().add(pdot2);

							Sphere pdot3 = new Sphere(1);
							pdot3.setMaterial(new PhongMaterial(edgecolor));
							pdot3.getTransforms().add(new Translate(points3dmatrix[a - 1][5].getX(),
									points3dmatrix[a - 1][5].getY(), points3dmatrix[a - 1][5].getZ()));
							group.getChildren().add(pdot3);

						}
					}
				}
			}
		}

		// Add lines
		if (false) {
			// Lines for barycentric subdivision of chambers:

			// for (int a = 1; a <= fDomain.size(); a++) {
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getVertex3D(0, a), fDomain.getEdgeCenter3D(1, a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getEdgeCenter3D(1, a), fDomain.getVertex3D(2, a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			//
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getVertex3D(2, a), fDomain.getEdgeCenter3D(0, a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getEdgeCenter3D(0, a), fDomain.getVertex3D(1, a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			//
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getVertex3D(0, a), fDomain.getChamberCenter3D(a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getChamberCenter3D(a), fDomain.getEdgeCenter3D(0, a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			//
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getVertex3D(1, a), fDomain.getChamberCenter3D(a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getChamberCenter3D(a), fDomain.getEdgeCenter3D(1, a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			//
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getVertex3D(2, a), fDomain.getChamberCenter3D(a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f));
			// group.getChildren().add(Cylinderline.createConnection(fDomain.
			// getChamberCenter3D(a), fDomain.getEdgeCenter3D(2, a),
			// Color.WHITE.deriveColor(0, 1, 1, 0.4), 0.5f)); }

			double width = 0;
			if (fDomain.getGeometry() == Geometry.Hyperbolic) {
				Point3D refPoint = fDomain.getChamberCenter3D(1).multiply(0.01);
				Point3D origin = new Point3D(0, 0, 1);
				double w = 0.01;
				double h = (1 + w * w) / (1 - w * w);
				// Length of translation
				double t = Tools.distance(fDomain, refPoint, origin);
				// Affine translation:
				Affine translateT = new Affine(Math.cosh(t), Math.sinh(t), 0, Math.sinh(t), Math.cosh(t), 0); // Translation
																												// along
																												// x-axis
				Point2D x = translateT.transform(0, 1);
				Point2D y = translateT.transform((1 + h) * w, h);

				width = 100 * (y.getX() / (1 + y.getY()) - x.getX() / (1 + x.getY()));
			} else if (fDomain.getGeometry() == Geometry.Euclidean) {
				width = 1;
			} else if (fDomain.getGeometry() == Geometry.Spherical) {
				width = 0.5;
			}

			// Edges of Tiling:
			Point3D v0, e2, v1;
			int m = fDomain.size();
			BitSet visited = new BitSet(m); //
			int a = 1;
			// Fallunterscheidung needs some serious refactoring
			if (geom == Geometry.Euclidean) {
				while (a <= m) {
					if (!visited.get(a)) {

						v0 = fDomain.getVertex3D(0, a);
						v1 = fDomain.getVertex3D(1, a);
						e2 = fDomain.getEdgeCenter3D(2, a);

						group.getChildren().add(Cylinderline.createConnection(v0, e2, Color.BLACK, width));
						group.getChildren().add(Cylinderline.createConnection(e2, v1, Color.BLACK, width));
						visited.set(a);
						visited.set(dsymbol.getS2(a));
					}
					a++;
				}
			} else if (geom == Geometry.Hyperbolic) {
				while (a <= m) {
					if (!visited.get(a)) {
						v0 = fDomain.getVertex3D(0, a);
						e2 = fDomain.getEdgeCenter3D(2, a);
						v1 = fDomain.getVertex3D(1, a);

						Point3D[] linePoints = new Point3D[9];
						linePoints[0] = v0;
						linePoints[4] = e2;
						linePoints[8] = v1;
						for (int i = 1; i < 4; i++) {
							linePoints[i] = Tools.interpolateHyperbolicPoints(v0, e2, i / 8d);
						}
						for (int i = 5; i < 8; i++) {
							linePoints[i] = Tools.interpolateHyperbolicPoints(e2, v1, i / 8d);
						}
						for (int i = 0; i < 8; i++) {
							group.getChildren().add(Cylinderline.createConnection(linePoints[i], linePoints[i + 1],
									Color.BLACK, width));
						}
						visited.set(dsymbol.getS2(a));
					}
					a++;
				}

			} else {
				// performanceProbleme durch Adden, also Magic numbers
				while (a <= m) {
					if (!visited.get(a)) {
						v0 = fDomain.getVertex3D(0, a);
						e2 = fDomain.getEdgeCenter3D(2, a);
						v1 = fDomain.getVertex3D(1, a);

						Point3D[] linePoints = new Point3D[33];
						linePoints[0] = v0;
						linePoints[16] = e2;
						linePoints[32] = v1;
						for (int i = 1; i < 16; i++) {
							linePoints[i] = Tools.interpolateSpherePoints(v0, e2, i / 32.0);
						}
						for (int i = 17; i < 32; i++) {
							linePoints[i] = Tools.interpolateSpherePoints(e2, v1, i / 32.0);
						}
						for (int j = 0; j < 32; j++) {
							group.getChildren().add(Cylinderline.createConnection(linePoints[j], linePoints[j + 1],
									Color.BLACK, width));
						}
						visited.set(dsymbol.getS2(a));
					}
					a++;
				}
			}

		}

		// add numbers:
		if (false) {
			for (int a = 1; a <= fDomain.size(); a++) {
				final Point3D apt = fDomain.getChamberCenter3D(a);
				Text label = new Text("" + a);
				label.setFont(Font.font(8));
				label.getTransforms().add(new Translate(apt.getX() - 4, apt.getY() + 4, apt.getZ()));

				label.setFill(Color.BLACK.deriveColor(0, 1, 1, 0.4));
				group.getChildren().add(label);
			}
		}

		// add some points to debug transforms:

		if (false) {
			for (int i = 0; i < 3; i++) {
				final Point3D a = fDomain.getVertex3D(i, 16);
				final Sphere sphere = new Sphere(2);
				switch (i) {
				case 0:
					sphere.setMaterial(new PhongMaterial(Color.GREEN));
					break;
				case 1:
					sphere.setMaterial(new PhongMaterial(Color.YELLOW));
					break;
				case 2:
					sphere.setMaterial(new PhongMaterial(Color.RED));
					break;
				}
				sphere.getTransforms().add(new Translate(a.getX(), a.getY(), a.getZ()));
				group.getChildren().add(sphere);
			}

			final Transform transform = Tiling.getTransform(fDomain.getGeometry(), fDomain.getVertex3D(0, 16),
					fDomain.getVertex3D(1, 16), fDomain.getVertex3D(0, 19), fDomain.getVertex3D(1, 19), true);

			for (int i = 0; i < 3; i++) {
				final Point3D a = fDomain.getVertex3D(i, 16);
				final Sphere sphere = new Sphere(2);
				sphere.getTransforms().addAll(transform, new Translate(a.getX(), a.getY(), a.getZ()));

				switch (i) {
				case 0:
					sphere.setMaterial(new PhongMaterial(Color.LIGHTGREEN));
					break;
				case 1:
					sphere.setMaterial(new PhongMaterial(Color.LIGHTYELLOW));
					break;
				case 2:
					sphere.setMaterial(new PhongMaterial(Color.PINK));
					break;

				}
				group.getChildren().add(sphere);
			}
		}
		return group;
	}

	// additional functions

	public static int[] invertOrientation(int[] arr) {
		int[] invArr = Arrays.copyOf(arr, arr.length);
		for (int i = 0; i < invArr.length / 6; i++) {
			int save = invArr[i * 6 + 2];
			invArr[i * 6 + 2] = invArr[i * 6 + 4];
			invArr[i * 6 + 4] = save;
		}

		return invArr;
	}

	public static double computeWindingNumber(Point3D a0, Point3D a1, Point3D a2) {
		return (a1.getX() - a0.getX()) * (a1.getY() + a0.getY()) + (a2.getX() - a1.getX()) * (a2.getY() + a1.getY())
				+ (a0.getX() - a2.getX()) * (a0.getY() + a2.getY());
	}

	public static double getLong(Point3D a) {
		return Math.atan2(a.getY(), a.getX());
	}

	public static double getLang(Point3D a, double radius) {
		return Math.asin(a.getZ() / radius);
	}

	public static double sphericalDistance(Point3D a, Point3D b, double radius) {
		return radius * Math.acos(a.normalize().dotProduct(b.normalize()));
	}

	public static double[][] sort2DArrayBasedOnSecondColumn(double[][] array) {
		Arrays.sort(array, Comparator.comparingDouble(arr -> arr[1]));
		return array;
	}

	public static Point3D getSphericalNormal(Point3D point) {
		Point3D p = new Point3D(2 * point.getX(), 2 * point.getY(), 2 * point.getZ()).normalize();
		return p;
	}

	public static Point3D getHyperbolicNormal(Point3D point) {
		Point3D p = new Point3D(2 * point.getX(), 2 * point.getY(), -2 * point.getZ()).normalize();
		return p;
	}

}
