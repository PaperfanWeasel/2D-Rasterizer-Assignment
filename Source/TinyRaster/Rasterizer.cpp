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
#include <stdlib.h>

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

	return true;
}

void Rasterizer::WriteRGBAToFramebuffer(int x, int y, const Colour4 & colour)
{
	if (x >= 0 && x < mWidth && y >= 0 && y < mHeight) 
	{
		PixelRGBA *pixel = mFramebuffer->GetBuffer();
	
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
	Vector2 pt1 = v1.position;
	Vector2 pt2 = v2.position;

	Colour4 colour1 = v1.colour;
	Colour4 colour2 = v2.colour;

	bool swap_vertices = pt1[0] > pt2[0]; // Case 1: set to true when x1 > x2

	if (swap_vertices)
	{
		pt1 = v2.position;
		pt2 = v1.position;
	}

	bool negative_slope = (pt2[1] - pt1[1]) < 0; // Case 2: set to true if y2-y1<0

	float dx = pt2[0] - pt1[0];
	float dy = pt2[1] - pt1[1];

	int epsilon = 0;

	int x = pt1[0];
	int y = pt1[1];
	int ex = pt2[0];
	int ey = pt2[1];

	bool swap_xy = abs(dx) < abs(dy); // Case 3: |dx| < |dy|

	if (swap_xy)
	{
		x = pt1[1]; 
		y = pt1[0]; 
		ex = pt2[1];
		ey = pt2[0];

		dy = pt2[0] - pt1[0];
		dx = pt2[1] - pt1[1];
	}

	if (negative_slope)
	{
		y = -y; // reflect the line about x axis
		dy = -dy; // also negate dy
	}
	// Special cases for octants 3 and 7
	if (negative_slope && swap_xy)
	{
		if (swap_vertices)
		{
			// Octant 3 Case: swap first then reflect
			pt1 = v1.position;
			pt2 = v2.position;

			x = pt1[1];
			y = pt1[0];
			ex = pt2[1];
			ey = pt2[0];
			dy = pt2[0] - pt1[0];
			dx = pt2[1] - pt1[1];

			y = -y;
			dy = -dy;
		}
		else 
		{
			// Octant 7 Case: reflect first then swap
			pt1[1] = -pt1[1];
			pt2[1] = -pt2[1];
			dy = -dy;

			x = pt1[1];
			y = pt1[0];
			ex = pt2[1];
			ey = pt2[0];
			dx = pt2[1] - pt1[1];
			dy = pt2[0] - pt1[0];
		}
	}
	while (x <= ex)
	{
		Vector2 temp(x, y);
		Colour4 colour = v1.colour;
		if (swap_xy)
		{
			temp[0] = y;
			temp[1] = x;
		}
		if (negative_slope)
		{
			temp[1] = -temp[1];
		}
		if (negative_slope && swap_xy)
		{
			if (swap_vertices)
			{
				temp[0] = -y;
				temp[1] = x;
			}
			else
			{
				temp[0] = y;
				temp[1] = -x;
			}
		}
		// Linear Interpolation
		if (mFillMode == INTERPOLATED_FILLED)
		{
			// t = |P0 - P| / |P1 - P0|
			Vector2 position = (v1.position - temp);
			Vector2 length = (v1.position - v2.position);
			float normPosition = position.Norm();
			float normLength = length.Norm();
			float t = normPosition / normLength;

			// C = t * C2 + ( 1 - t ) * C1
			colour = colour2*t + colour1*(1.0 - t);
		}
		SetFGColour(colour);
		// Task 3: Handling Thickness
		if (thickness > 1)
		{
			int j = thickness / 2;
			for (int i = 1; i < j; i++)
			{
				if (swap_xy) // check if line is horizontal or vertical
				{
					DrawPoint2D(Vector2(temp[0] - i, temp[1])); // left
					DrawPoint2D(Vector2(temp[0] + i, temp[1])); // right
				}
				else
				{
					DrawPoint2D(Vector2(temp[0], temp[1] - i)); // below
					DrawPoint2D(Vector2(temp[0], temp[1] + i)); // above
				}
			}
		}
		DrawPoint2D(temp);

		epsilon += dy;

		if ((epsilon << 1) >= dx)
		{
			y++;
			
			epsilon -= dx;
		}

		x++;
	}
}

void Rasterizer::DrawUnfilledPolygon2D(const Vertex2d * vertices, int count)
{
	for (int i = 0; i < (count-1); i++)
	{
		DrawLine2D(vertices[i], vertices[i + 1]);
	}
	DrawLine2D(vertices[count-1], vertices[0]);
}

bool sortCheckFunction(ScanlineLUTItem i, ScanlineLUTItem j) { return (i.pos_x<j.pos_x); }

void Rasterizer::ScanlineFillPolygon2D(const Vertex2d * vertices, int count)
{
	ClearScanlineLUT(); // Clear look up table

	for (int i = 0; i < count; i++)
	{
		Vector2 pt1 = vertices[i].position;
		Vector2 pt2 = vertices[(i + 1) % count].position;
		Colour4 colour = vertices[i].colour;
		float x1 = pt1[0];
		float y1 = pt1[1];
		float x2 = pt2[0];
		float y2 = pt2[1];
		
		// populating look up table
		for (int y = 0; y < mHeight - 1; y++)
		{
			// check to see if scanline intersects with polygon
			if ((y >= y1 && y <= y2) || (y <= y1 && y >= y2))
			{
				if (y1 != y2) // Horizontal Line detection
				{
					float m = ((x1 - x2) / (y1 - y2)); // Gradient
					float x = x1 + (y - y1) * m;
					// Add item to the look up table
					ScanlineLUTItem newitem = { colour , x };
					mScanlineLUT[y].push_back(newitem);
				}
			}
		}
	}
	// Filling Scanlines
	for (int i = 0; i < mHeight - 1; i++)
	{
		std::sort(mScanlineLUT[i].begin(), mScanlineLUT[i].end(), sortCheckFunction); // Sorting using sortCheckFunction
		if (!mScanlineLUT[i].empty())
		{
			if (mScanlineLUT[i].size() % 2 == 0) // check if even
			{
				for (int j = 0; j < mScanlineLUT[i].size(); j+=2 )
				{
					// two adjacent points in the LUT on the i-th scanline
					Vector2 pt1(mScanlineLUT[i][j].pos_x, i);
					Vector2 pt2(mScanlineLUT[i][j + 1].pos_x, i);
					
					Vertex2d v1{ mScanlineLUT[i][j].colour, pt1 };
					Vertex2d v2{ mScanlineLUT[i][j + 1].colour, pt2 };

					DrawLine2D(v1, v2);
				}
			}
		}
	}
	//Ex 2.3 Extend Rasterizer::ScanlineFillPolygon2D method so that it is capable of alpha blending, i.e. draw translucent polygons.
	//Note: The variable mBlendMode indicates if the blend mode is set to alpha blending.
	//To do alpha blending during filling, the new colour of a point should be combined with the existing colour in the framebuffer using the alpha value.
	//Use Test 6 (Press F6) to test your solution
}

void Rasterizer::ScanlineInterpolatedFillPolygon2D(const Vertex2d * vertices, int count)
{
	ClearScanlineLUT(); // Clear look up table
	if (mFillMode == INTERPOLATED_FILLED)
	{
		for (int i = 0; i < count; i++)
		{
			Vector2 pt1 = vertices[i].position;
			Vector2 pt2 = vertices[(i + 1) % count].position;
			Colour4 colour1 = vertices[i].colour;
			Colour4 colour2 = vertices[(i + 1) % count].colour;
			float x1 = pt1[0];
			float y1 = pt1[1];
			float x2 = pt2[0];
			float y2 = pt2[1];

			// populating look up table
			for (int y = 0; y < mHeight - 1; y++)
			{
				// check to see if scanline intersects with polygon
				if ((y >= y1 && y <= y2) || (y <= y1 && y >= y2))
				{
					if (y1 != y2) // Horizontal Line detection
					{
						float m = ((x1 - x2) / (y1 - y2)); // Gradient
						float x = x1 + (y - y1) * m;
						Vector2 temp(x, y);
						// Colour Interpolation
						Vector2 position = (pt1 - temp);
						Vector2 length = (pt1 - pt2);
						float normPosition = position.Norm();
						float normLength = length.Norm();
						float t = normPosition / normLength;
						Colour4 colour = colour2*t + colour1*(1.0 - t);
						// Add item to the look up table
						ScanlineLUTItem newitem = { colour , x };
						mScanlineLUT[y].push_back(newitem);
					}
				}
			}
		}
		// Filling Scanlines
		for (int i = 0; i < mHeight - 1; i++)
		{
			std::sort(mScanlineLUT[i].begin(), mScanlineLUT[i].end(), sortCheckFunction); // Sorting using sortCheckFunction
			if (!mScanlineLUT[i].empty())
			{
				if (mScanlineLUT[i].size() % 2 == 0)
				{
					for (int j = 0; j < mScanlineLUT[i].size(); j += 2)
					{
						// two adjacent points in the LUT on the i-th scanline
						Vector2 pt1(mScanlineLUT[i][j].pos_x, i);
						Vector2 pt2(mScanlineLUT[i][j + 1].pos_x, i);

						Vertex2d v1{ mScanlineLUT[i][j].colour, pt1 };
						Vertex2d v2{ mScanlineLUT[i][j + 1].colour, pt2 };

						DrawLine2D(v1, v2);
					}
				}
			}
		}
	}
}

void Rasterizer::DrawCircle2D(const Circle2D & inCircle, bool filled)
{
	int r = inCircle.radius;
	Vector2 c = inCircle.centre;
	Colour4 colour = inCircle.colour;
	SetFGColour(colour);

	int cx, cy, r2, x, y;
	
	cx = c[0];
	cy = c[1];
	r2 = r * r;
	// Drawing Circles using 8 way symmetry
	if (filled)
	{
		DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx, cy + r) } , 5);
		DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx, cy - r) } , 5);
		DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx + r, cy) } , 5);
		DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx - r, cy) } , 5);
		x = 0;
		y = (int)(sqrt(r2 - 1) + 0.5);
		while (x <= y)
		{
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx + x, cy + y) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx + x, cy - y) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx - x, cy + y) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx - x, cy - y) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx + y, cy + x) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx + y, cy - x) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx - y, cy + x) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx - y, cy - x) } , 5);
			x++;
			y = (int)(sqrt(r2 - x*x) + 0.5);
		}
		if (x == y)
		{
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx + x, cy + y) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx + x, cy - y) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx - x, cy + y) } , 5);
			DrawLine2D(Vertex2d{ colour, c }, Vertex2d{ colour, Vector2(cx - x, cy - y) } , 5);
		}
	}
	else // Draw unfilled
	{
		DrawPoint2D(Vector2(cx, cy + r));
		DrawPoint2D(Vector2(cx, cy - r));
		DrawPoint2D(Vector2(cx + r, cy));
		DrawPoint2D(Vector2(cx - r, cy));
		x = 0;
		y = (int)(sqrt(r2 - 1) + 0.5);
		while (x < y)
		{
			DrawPoint2D(Vector2(cx + x, cy + y));
			DrawPoint2D(Vector2(cx + x, cy - y));
			DrawPoint2D(Vector2(cx - x, cy + y));
			DrawPoint2D(Vector2(cx - x, cy - y));
			DrawPoint2D(Vector2(cx + y, cy + x));
			DrawPoint2D(Vector2(cx + y, cy - x));
			DrawPoint2D(Vector2(cx - y, cy + x));
			DrawPoint2D(Vector2(cx - y, cy - x));
			x++;
			y = (int)(sqrt(r2 - x*x) + 0.5);
		}
		if (x == y)
		{
			DrawPoint2D(Vector2(cx + x, cy + y));
			DrawPoint2D(Vector2(cx + x, cy - y));
			DrawPoint2D(Vector2(cx - x, cy + y));
			DrawPoint2D(Vector2(cx - x, cy - y));
		}
	}
}

Framebuffer *Rasterizer::GetFrameBuffer() const
{
	return mFramebuffer;
}
