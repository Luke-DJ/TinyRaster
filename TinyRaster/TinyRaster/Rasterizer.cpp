/*---------------------------------------------------------------------
*
* Copyright © 2016  Minsi Chen
* E-mail: m.chen@derby.ac.uk
*
* The source is written for the Graphics I and II modules. You are free
* to use and extend the functionality. The code provided here is functional
* however the author does not guarantee its performance.
---------------------------------------------------------------------*/
#include <algorithm>
#include <math.h>

#include "Rasterizer.h"

Rasterizer::Rasterizer(void)
{
	mFramebuffer = NULL;
	mScanlineLUT = NULL;
}

void Rasterizer::ClearScanlineLUT()
{
	Scanline *pScanline = mScanlineLUT;

	for (int y = 0; y < mHeight; y++)
	{
		(pScanline + y)->clear();
		(pScanline + y)->shrink_to_fit();
	}
}

unsigned int Rasterizer::ComputeOutCode(const Vector2 & p, const ClipRect& clipRect)
{
	unsigned int CENTRE = 0x0;
	unsigned int LEFT = 0x1;
	unsigned int RIGHT = 0x1 << 1;
	unsigned int BOTTOM = 0x1 << 2;
	unsigned int TOP = 0x1 << 3;
	unsigned int outcode = CENTRE;
	
	if (p[0] < clipRect.left)
		outcode |= LEFT;
	else if (p[0] >= clipRect.right)
		outcode |= RIGHT;

	if (p[1] < clipRect.bottom)
		outcode |= BOTTOM;
	else if (p[1] >= clipRect.top)
		outcode |= TOP;

	return outcode;
}

bool Rasterizer::ClipLine(const Vertex2d & v1, const Vertex2d & v2, const ClipRect& clipRect, Vector2 & outP1, Vector2 & outP2)
{
	//TODO: EXTRA This is not directly prescribed as an assignment exercise. 
	//However, if you want to create an efficient and robust rasteriser, clipping is a usefull addition.
	//The following code is the starting point of the Cohen-Sutherland clipping algorithm.
	//If you complete its implementation, you can test it by calling prior to calling any DrawLine2D .

	const Vector2 p1 = v1.position;
	const Vector2 p2 = v2.position;
	unsigned int outcode1 = ComputeOutCode(p1, clipRect);
	unsigned int outcode2 = ComputeOutCode(p2, clipRect);

	outP1 = p1;
	outP2 = p2;
	
	bool draw = false;

	return true;
}

void Rasterizer::WriteRGBAToFramebuffer(int x, int y, const Colour4 & colour)
{
	//Ensuring that pixels outside of the framebuffer region aren't added to the framebuffer
	if (x < 0 || x >= mWidth || y < 0 || y >= mHeight)
	{
		return;
	}

	PixelRGBA *pixel = mFramebuffer->GetBuffer();
	
	/*If mBlendMode is set to ALPHA_BLEND, combine the new colour (colour) with the old colour (pixel[y*mWidth + x])
	Otherwise, set the current colour to the new colour*/
	if (mBlendMode == ALPHA_BLEND)
	{
		pixel[y*mWidth + x] = (colour * colour[3]) + (pixel[y*mWidth + x] * (1 - colour[3]));
	}
	else
	{
		pixel[y*mWidth + x] = colour;
	}
}

Rasterizer::Rasterizer(int width, int height)
{
	//Initialise the rasterizer to its initial state
	mFramebuffer = new Framebuffer(width, height);
	mScanlineLUT = new Scanline[height];
	mWidth = width;
	mHeight = height;

	mBGColour.SetVector(0.0, 0.0, 0.0, 1.0);	//default bg colour is black
	mFGColour.SetVector(1.0, 1.0, 1.0, 1.0);    //default fg colour is white

	mGeometryMode = LINE;
	mFillMode = UNFILLED;
	mBlendMode = NO_BLEND;

	SetClipRectangle(0, mWidth, 0, mHeight);
}

Rasterizer::~Rasterizer()
{
	delete mFramebuffer;
	delete[] mScanlineLUT;
}

void Rasterizer::Clear(const Colour4& colour)
{
	PixelRGBA *pixel = mFramebuffer->GetBuffer();

	SetBGColour(colour);

	int size = mWidth*mHeight;
	
	for(int i = 0; i < size; i++)
	{
		//fill all pixels in the framebuffer with background colour
		*(pixel + i) = mBGColour;
	}
}

void Rasterizer::DrawPoint2D(const Vector2& pt, int size)
{
	int x = pt[0];
	int y = pt[1];

	WriteRGBAToFramebuffer(x, y, mFGColour);
}

void Rasterizer::DrawLine2D(const Vertex2d & v1, const Vertex2d & v2, int thickness)
{
	//The following code is basic Bresenham's line drawing algorithm.
	//The current implementation is only capable of rasterise a line in the first octant, where dy < dx and dy/dx >= 0
	//See if you want to read ahead https://www.cs.helsinki.fi/group/goa/mallinnus/lines/bresenh.html
	
	//TODO:
	//Ex 1.1 Complete the implementation of Rasterizer::DrawLine2D method. 
	//This method currently consists of a partially implemented Bresenham algorithm.
	//You must extend its implementation so that it is capable of drawing 2D lines with arbitrary gradient(slope).
	//Use Test 1 (Press F1) to test your implementation
	
	//Ex 1.2 Extend the implementation of Rasterizer::DrawLine2D so that it is capable of drawing lines based on a given thickness.
	//Note: The thickness is passed as an argument int thickness.
	//Use Test 2 (Press F2) to test your implementation

	//Ex 1.3 Extend the implementation of Rasterizer::DrawLine2D so that it is capable of interpolating colour across a line when each end-point has different colours.
	//Note: The variable mFillMode indicates if the fill mode is set to INTERPOLATED_FILL. 
	//The colour of each point should be linearly interpolated using the colours of v1 and v2.
	//Use Test 2 (Press F2) to test your implementation

	//Refactored from 'drawline_Bresenham' in 'ASCiiDraw'

	Vector2 pt1 = v1.position;
	Vector2 pt2 = v2.position;

	bool swap_x = pt2[0] < pt1[0];

	int dx = swap_x ? pt1[0] - pt2[0] : pt2[0] - pt1[0];
	int dy = swap_x ? pt1[1] - pt2[1] : pt2[1] - pt1[1];

	int reflect = dy < 0 ? -1 : 1;
	bool swap_xy = dy*reflect > dx;
	int epsilon = 0;

	int sx = swap_xy ? reflect < 0 ? swap_x ? pt1[1] : pt2[1] : swap_x ? pt2[1] : pt1[1] : swap_x ? pt2[0] : pt1[0];
	int y = swap_xy ? reflect < 0 ? swap_x ? pt1[0] : pt2[0] : swap_x ? pt2[0] : pt1[0] : swap_x ? pt2[1] : pt1[1];
	int ex = swap_xy ? reflect < 0 ? swap_x ? pt2[1] : pt1[1] : swap_x ? pt1[1] : pt2[1] : swap_x ? pt1[0] : pt2[0];

	int x = sx;
	y *= reflect;

	while (x <= ex)
	{
		Vector2 temp(swap_xy ? y*reflect : x, swap_xy ? x : y*reflect);

		//Adding colour
		Colour4 colour;
		float t;
		if (swap_xy)
		{
			t = abs((pt1[1] - temp[1]) / (pt2[1] - pt1[1]));
		}
		else
		{
			t = abs((pt1[0] - temp[0]) / (pt2[0] - pt1[0]));
		}
		colour = (v2.colour * t) + (v1.colour * (1 - t));
		SetFGColour(colour);

		//Drawing the current point
		DrawPoint2D(Vector2(temp[0], temp[1]));

		//Distributing 'thickness' evenly on each side of the current point
		if (swap_xy)
		{
			for (int i = 1; i <= ceil((thickness - 1) / 2.0); i++)
			{
				DrawPoint2D(Vector2(temp[0] + i, temp[1]));
			}
			for (int i = 1; i <= floor((thickness - 1) / 2.0); i++)
			{
				DrawPoint2D(Vector2(temp[0] - i, temp[1]));
			}
		}
		else
		{
			for (int i = 1; i <= ceil((thickness - 1) / 2.0); i++)
			{
				DrawPoint2D(Vector2(temp[0], temp[1] + i));
			}
			for (int i = 1; i <= floor((thickness - 1) / 2.0); i++)
			{
				DrawPoint2D(Vector2(temp[0], temp[1] - i));
			}
		}

		epsilon += swap_xy ? dx : dy * reflect;

		if ((epsilon << 1) >= (swap_xy ? dy * reflect : dx))
		{
			y++;

			epsilon -= swap_xy ? dy * reflect : dx;
		}
		x++;
	}
}

void Rasterizer::DrawUnfilledPolygon2D(const Vertex2d * vertices, int count)
{
	//TODO:
	//Ex 2.1 Implement the Rasterizer::DrawUnfilledPolygon2D method so that it is capable of drawing an unfilled polygon, i.e. only the edges of a polygon are rasterised. 
	//Please note, in order to complete this exercise, you must first complete Ex1.1 since DrawLine2D method is reusable here.
	//Note: The edges of a given polygon can be found by conntecting two adjacent vertices in the vertices array.
	//Use Test 3 (Press F3) to test your solution.

	//Drawing a line between each pair of vertices in the vertices array
	for (int i = 0; i < count; i++)
	{
		if (i + 1 == count)
		{
			DrawLine2D(vertices[i], vertices[0]);
		}
		else
		{
			DrawLine2D(vertices[i], vertices[i + 1]);
		}
	}
}

void Rasterizer::ScanlineFillPolygon2D(const Vertex2d * vertices, int count)
{
	//TODO:
	//Ex 2.2 Implement the Rasterizer::ScanlineFillPolygon2D method method so that it is capable of drawing a solidly filled polygon.
	//Note: You can implement floodfill for this exercise however scanline fill is considered a more efficient and robust solution.
	//		You should be able to reuse DrawUnfilledPolygon2D here.
	//
	//Use Test 4 (Press F4) to test your solution, this is a simple test case as all polygons are convex.
	//Use Test 5 (Press F5) to test your solution, this is a complex test case with one non-convex polygon.

	//Ex 2.3 Extend Rasterizer::ScanlineFillPolygon2D method so that it is capable of alpha blending, i.e. draw translucent polygons.
	//Note: The variable mBlendMode indicates if the blend mode is set to alpha blending.
	//To do alpha blending during filling, the new colour of a point should be combined with the existing colour in the framebuffer using the alpha value.
	//Use Test 6 (Press F6) to test your solution

	ClearScanlineLUT();

	/*Iterating over all points in the vertices array, for each scanline
	If one point is above/below, and the other is below/above, in a pair of points in the vertices array,
	then calculate the x position of the intersect and add it, along with the current point's colour/calculated colour, to the current y position of the LUT*/
	for (int i = 0; i < mHeight; i++)
	{
		for (int j = 0; j < count; j++)
		{
			int y = i;
			int x1 = vertices[j].position[0];
			int y1 = vertices[j].position[1];
			int x2;
			int y2;
			if (j + 1 == count)
			{
				x2 = vertices[0].position[0];
				y2 = vertices[0].position[1];
			}
			else
			{
				x2 = vertices[j + 1].position[0];
				y2 = vertices[j + 1].position[1];
			}

			if ((y >= y1 && y <= y2) || (y <= y1 && y >= y2))
			{
				int x;
				if (y1 - y2 != 0)
				{
					x = x1 + (y - y1) * (x1 - x2) / (y1 - y2);

					ScanlineLUTItem item;

					if (mFillMode == INTERPOLATED_FILLED)
					{
						float a = y1 - y;
						float b = y2 - y1;
						float t = abs(a / b);

						Colour4 colour;
						if (j + 1 == count)
						{
							colour = (vertices[0].colour * t) + (vertices[j].colour * (1 - t));
						}
						else
						{
							colour = (vertices[j + 1].colour * t) + (vertices[j].colour * (1 - t));
						}

						item = { colour, x };
					}
					else
					{
						item = { vertices[j].colour, x };
					}

					mScanlineLUT[y].push_back(item);
				}
			}
		}
	}
	//Removing duplicates from, and sorting, each scanline
	for (int i = 0; i < mHeight; i++)
	{
		for (int j = 0; j + 1 < mScanlineLUT[i].size(); j++)
		{
			if (mScanlineLUT[i][j].pos_x == mScanlineLUT[i][j + 1].pos_x)
			{
				mScanlineLUT[i].erase(begin(mScanlineLUT[i]) + j);
				j--;
			}
		}
		sort(begin(mScanlineLUT[i]), end(mScanlineLUT[i]), [](ScanlineLUTItem a, ScanlineLUTItem b) {return a.pos_x < b.pos_x; });
	}

	/*Creating a vector of vectors, that contain Vector2s, for the vertices that the filling algorithm should treat uniquely
	This is a 2d vector rather than a standard vector, so that when checking to see if a point exists in ignoreVertices,
	it's compared to points in ignoreVertices on the current Scanline, rather than every point in ignoreVertices*/
	std::vector<std::vector<Vector2>> ignoreVertices = {};
	for (int i = 0; i <= mHeight; i++)
	{
		ignoreVertices.push_back(std::vector<Vector2>());
	}
	/*Iterating over all points in the vertices array.
	If the previous point is above and the next point is above, or the previous point is below and the next point is below, add the current point to ignoreVertices*/
	for (int i = 0; i < count; i++)
	{
		if (vertices[i].position[1] >= 0 && vertices[i].position[1] <= mHeight)
		{
			if (i == 0)
			{
				if ((vertices[i].position[1] > vertices[count - 1].position[1] && vertices[i].position[1] > vertices[i + 1].position[1]) ||
					(vertices[i].position[1] < vertices[count - 1].position[1] && vertices[i].position[1] < vertices[i + 1].position[1]))
				{
					ignoreVertices[vertices[i].position[1]].push_back(vertices[i].position);
				}
			}
			else if (i + 1 == count)
			{
				if ((vertices[i].position[1] > vertices[i - 1].position[1] && vertices[i].position[1] > vertices[0].position[1]) ||
					(vertices[i].position[1] < vertices[i - 1].position[1] && vertices[i].position[1] < vertices[0].position[1]))
				{
					ignoreVertices[vertices[i].position[1]].push_back(vertices[i].position);
				}
			}
			else
			{
				if ((vertices[i].position[1] > vertices[i - 1].position[1] && vertices[i].position[1] > vertices[i + 1].position[1]) ||
					(vertices[i].position[1] < vertices[i - 1].position[1] && vertices[i].position[1] < vertices[i + 1].position[1]))
				{
					ignoreVertices[vertices[i].position[1]].push_back(vertices[i].position);
				}
			}
		}
	}

	//Drawing lines between pairs of ScanlineLUTItems in each Scanline in the LUT
	for (int i = 0; i < mHeight; i++)
	{
		if (mScanlineLUT[i].size() > 1)
		{
			bool oddEvenBool = true;
			for (int j = 0; j + 1 < mScanlineLUT[i].size(); j++)
			{
				bool ignoreVertex = false;
				//Setting ignoreVertex to true if the current Vertex2d's position is in ignoreVertices
				for (int k = 0; k < ignoreVertices[i].size(); k++)
				{
					if (ignoreVertices[i][k][0] == mScanlineLUT[i][j].pos_x)
					{
						ignoreVertex = true;
						break;
					}
				}

				//Drawing a point at j + 1 if the last item in mScanlineLUT[i] is in ignoreVertices
				if (j + 2 == mScanlineLUT[i].size())
				{
					for (int k = 0; k < ignoreVertices[i].size(); k++)
					{
						if (ignoreVertices[i][k][0] == mScanlineLUT[i][j + 1].pos_x)
						{
							SetFGColour(mScanlineLUT[i][j + 1].colour);
							DrawPoint2D(Vector2(mScanlineLUT[i][j + 1].pos_x, i));
						}
					}
				}

				//Ensuring the position of items in ignoreVertices are drawn if they're supposed to be
				if (ignoreVertex && oddEvenBool)
				{
					SetFGColour(mScanlineLUT[i][j].colour);
					DrawPoint2D(Vector2(mScanlineLUT[i][j].pos_x, i));
				}
				/*if ignoreVertex is true and oddEvenBool is false, or ignoreVertex is false and oddEvenBool is true,
				draw a line between the current ScanlineLUTItem and the next one*/
				else if (ignoreVertex != oddEvenBool)
				{
					Vertex2d vertexOne = Vertex2d{ mScanlineLUT[i][j].colour, Vector2(mScanlineLUT[i][j].pos_x, i) };
					Vertex2d vertexTwo = Vertex2d{ mScanlineLUT[i][j + 1].colour, Vector2(mScanlineLUT[i][j + 1].pos_x, i) };
					DrawLine2D(vertexOne, vertexTwo);
					oddEvenBool = false;
				}
				else
				{
					oddEvenBool = true;
				}
			}
		}
		//If there's only one ScanlineLUTItem in the Scanline, draw at that ScanlineLUTItem's position
		else if (mScanlineLUT[i].size() == 1)
		{
			SetFGColour(mScanlineLUT[i][0].colour);
			DrawPoint2D(Vector2(mScanlineLUT[i][0].pos_x, i));
		}
	}
}

void Rasterizer::ScanlineInterpolatedFillPolygon2D(const Vertex2d * vertices, int count)
{
	//TODO:
	//Ex 2.4 Implement Rasterizer::ScanlineInterpolatedFillPolygon2D method so that it is capable of performing interpolated filling.
	//Note: mFillMode is set to INTERPOLATED_FILL
	//		This exercise will be more straightfoward if Ex 1.3 has been implemented in DrawLine2D
	//Use Test 7 to test your solution

	/*Since ScanlineFillPolygon2D already takes care of filling,
	the interpolation code has been added to ScanlineFillPolygon2D so that code isn't needlessly repeated here*/

	ScanlineFillPolygon2D(vertices, count);
}

void Rasterizer::DrawCircle2D(const Circle2D & inCircle, bool filled)
{
	//TODO:
	//Ex 2.5 Implement Rasterizer::DrawCircle2D method so that it can draw a filled circle.
	//Note: For a simple solution, you can first attempt to draw an unfilled circle in the same way as drawing an unfilled polygon.
	//Use Test 8 to test your solution

	//How many points there will be on the circle's circumference (Pi x diameter) - 1 pixel step
	const int nSegment = acos(-1) * (inCircle.radius * 2);

	//The vector for the Vertex2ds of the circle's vertices
	std::vector<Vertex2d> verticesVector = {};

	//2 Pi Radians (360 degrees) divided by the number of vertices to get the angle of each 'slice' of the circle
	float slice = (2 * acos(-1)) / nSegment;

	/*If the number of segments is a multiple of 8,
	generate an eighth of the vertices,
	and add a new vertex to verticesVector of each existing vertex in reverse order,
	with the gap between the y/x of the original vertex and the y/x of the original final generated vertex added to the x/y of the new vertex, in order to form a quarter of the circle*/
	if (nSegment % 8 == 0)
	{
		for (int i = 0; i < (nSegment / 8) + 1; i++)
		{
			Vertex2d vertex = Vertex2d{
				inCircle.colour,
				Vector2(
					(int(inCircle.radius * cos(i * slice)) + inCircle.centre[0]),
					(int(inCircle.radius * sin(i * slice)) + inCircle.centre[1]))
			};
			//Ensuring only the outmost pixels of each line are added to verticesVector if the circle is to be filled
			if (filled)
			{
				if (vertex.position[0] >= inCircle.centre[0] && int(inCircle.radius * sin((i - 1) * slice)) + inCircle.centre[1] < vertex.position[1])
				{
					verticesVector.push_back(vertex);
				}
			}
			else
			{
				verticesVector.push_back(vertex);
			}
		}
		int firstVertices = verticesVector.size() - 1;
		for (int i = 1; i <= firstVertices; i++)
		{
			Vertex2d vertex = Vertex2d{
				inCircle.colour,
				Vector2(
					verticesVector[firstVertices].position[0] + (verticesVector[firstVertices - i].position[1] - verticesVector[firstVertices].position[1]),
					verticesVector[firstVertices].position[1] + (verticesVector[firstVertices - i].position[0] - verticesVector[firstVertices].position[0]))
			};
			//Ensuring only the outmost pixels of each line are added to verticesVector if the circle is to be filled
			if (filled)
			{
				if (verticesVector[verticesVector.size() - 1].position[1] < vertex.position[1])
				{
					verticesVector.push_back(vertex);
				}
			}
			else
			{
				verticesVector.push_back(vertex);
			}
		}
	}
	/*If the number of segments is a multiple of 4,
	generate a quarter of the vertices (if vertices weren't already generated due to the number of segments being a multiple of eight (or therefore potentially an even greater power of two)),
	and add a new vertex to verticesVector of each existing vertex in reverse order, with its x flipped to the opposite side of the centre position, in order to form a half of the circle*/
	if (nSegment % 4 == 0)
	{
		if (nSegment % 8 != 0)
		{
			for (int i = 0; i < (nSegment / 4) + 1; i++)
			{
				Vertex2d vertex = Vertex2d{
					inCircle.colour,
					Vector2(
						int(inCircle.radius * cos(i * slice)) + inCircle.centre[0],
						int(inCircle.radius * sin(i * slice)) + inCircle.centre[1])
				};
				//Ensuring only the outmost pixels of each line are added to verticesVector if the circle is to be filled
				if (filled)
				{
					if (int(inCircle.radius * sin((i - 1) * slice)) + inCircle.centre[1] < vertex.position[1])
					{
						verticesVector.push_back(vertex);
					}
				}
				else
				{
					verticesVector.push_back(vertex);
				}
			}
		}
		int firstVertices = verticesVector.size() - 1;
		for (int i = 1; i <= firstVertices; i++)
		{
			Vertex2d vertex = Vertex2d{
				inCircle.colour,
				Vector2(
					inCircle.centre[0] - (verticesVector[firstVertices - i].position[0] - inCircle.centre[0]),
					verticesVector[firstVertices - i].position[1])
			};
			verticesVector.push_back(vertex);
		}
	}
	/*Generate half of the vertices (if vertices weren't already generated due to the number of segments being a multiple of four (or therefore potentially an even greater power of two)),
	and add a new vertex to verticesVector of each existing vertex in reverse order, with its y flipped to the opposite side of the centre position, in order to form the whole circle*/
	if (nSegment % 4 != 0)
	{
		int firstVertices = nSegment / 2;
		if (nSegment % 2 == 0)
		{
			firstVertices++;
		}
		for (int i = 0; i < firstVertices; i++)
		{
			Vertex2d vertex = Vertex2d{
				inCircle.colour,
				Vector2(
					int(inCircle.radius * cos(i * slice)) + inCircle.centre[0],
					int(inCircle.radius * sin(i * slice)) + inCircle.centre[1])
			};
			//Ensuring only the outmost pixels of each line are added to verticesVector if the circle is to be filled
			if (filled)
			{
				if ((vertex.position[0] >= inCircle.centre[0] && int(inCircle.radius * sin((i - 1) * slice) + inCircle.centre[1]) < vertex.position[1]) ||
					(vertex.position[0] < inCircle.centre[0] && int(inCircle.radius * sin((i + 1) * slice) + inCircle.centre[1]) < vertex.position[1]))
				{
					verticesVector.push_back(vertex);
				}
			}
			else
			{
				verticesVector.push_back(vertex);
			}
		}
	}
	int secondVertices = verticesVector.size() - 1;
	for (int i = 1; i < secondVertices; i++)
	{
		Vertex2d vertex = Vertex2d{
			inCircle.colour,
			Vector2(
				verticesVector[secondVertices - i].position[0],
				inCircle.centre[1] - (verticesVector[secondVertices - i].position[1] - inCircle.centre[1]))
		};
		verticesVector.push_back(vertex);
	}

	//The array for the Vertex2ds (dynamically assigned to the size of verticesVector)
	Vertex2d* verticesArray = new Vertex2d[verticesVector.size()];
	copy(begin(verticesVector), end(verticesVector), verticesArray);

	//If filled is true, use ScanlineFillPolygon2D to draw a filled-in shape
	//Otherwise, use DrawUnfilledPolygon2D to draw an un-filled shape
	if (filled)
	{
		ScanlineFillPolygon2D(verticesArray, verticesVector.size());
	}
	else
	{
		DrawUnfilledPolygon2D(verticesArray, verticesVector.size());
	}

	delete[] verticesArray;
}

Framebuffer *Rasterizer::GetFrameBuffer() const
{
	return mFramebuffer;
}
