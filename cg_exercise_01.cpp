#include <cglib/colors/exercise.h>
#include <cglib/colors/convert.h>
#include <cglib/colors/cmf.h>

#include <cglib/core/glheaders.h>
#include <cglib/core/glmstream.h>

#include <cglib/core/assert.h>
#include <iostream>

/*
 * Teil-A
 * Draw the given vertices directly as GL_TRIANGLES.
 * For each vertex, also set the corresponding color.
 */

void draw_triangles(
	std::vector<glm::vec3> const& vertices,
	std::vector<glm::vec3> const& colors)
{
	cg_assert(vertices.size() == colors.size());
	cg_assert(vertices.size() % 3 == 0);


	glBegin(GL_TRIANGLES);
	for (int i = 0; i < vertices.size(); i++)
	{
		glColor3f(colors[i][0], colors[i][1], colors[i][2]);
		glVertex3f(vertices[i][0], vertices[i][1], vertices[i][2]);
	}
	glEnd();
}

/*
 * Teil-B
 * Generate a regular grid of resolution N*N (2*N*N triangles) in the xy-plane (z=0).
 * Store the grid in vertex-index form.
 *
 * The vertices of the triangles should be in counter clock-wise order.
 *
 * The grid must fill exactly the square [0, 1]x[0, 1], and must
 * be generated as an Indexed Face Set (Shared Vertex representation).
 *
 * The first vertex should be at position (0,0,0) and the last
 * vertex at position (1,1,0)
 *
 * An example for N = 3:
 *
 *   ^
 *   |  ----------
 *   |  |\ |\ |\ |
 *   |  | \| \| \|
 *   |  ----------
 *   |  |\ |\ |\ |
 * y |  | \| \| \|
 *   |  ----------
 *   |  |\ |\ |\ |
 *   |  | \| \| \|
 *   |  ----------
 *   |
 *   |-------------->
 *          x
 *
 */
void generate_grid(
	std::uint32_t N,
	std::vector<glm::vec3>* vertices,
	std::vector<glm::uvec3>* indices)
{
	cg_assert(N >= 1);
	cg_assert(vertices);
	cg_assert(indices);

	vertices->clear();
	indices->clear();

	std::uint32_t num_grid_vertices = (N+1) * (N+1);
	float stride = 1.0 / N;

	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			float x = j * stride;
			float y = i * stride;
			vertices->push_back(glm::vec3(x, y, 0.0));
		}
	}

	//indices of traingles
	for (std::uint32_t index = 0; index <= num_grid_vertices - N - 2; index++)
	{
		// only left_down_trangiles
		if (index % (N+1) == 0)
		{
			indices->push_back(glm::uvec3(index, index + 1, index + 1 + N));
		}
		// only right_up trangiles
		else if ((index+1) % (N+1) == 0)
		{
			indices->push_back(glm::uvec3(index, index + N + 1,index + N));
		}
		// both left_down and right_up trangiles
		else
		{
			indices->push_back(glm::uvec3(index, index + 1, index + 1 + N));
			indices->push_back(glm::uvec3(index, index + N + 1, index + N));
		}
	}
}

/*
 * Draw the given vertices as indexed GL_TRIANGLES using glDrawElements.
 * For each vertex, also set the corresponding color.
 *
 * Don't forget to enable the correct client states. After drawing
 * the triangles, you need to disable the client states again.
 */
void draw_indexed_triangles(
	std::vector<glm::vec3>  const& vertices,
	std::vector<glm::vec3>  const& colors,
	std::vector<glm::uvec3> const& indices)
{
	cg_assert(vertices.size() == colors.size());

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, vertices.data());
	glColorPointer(3, GL_FLOAT, 0, colors.data());

	glDrawElements(GL_TRIANGLES, indices.size() * 3, GL_UNSIGNED_INT, indices.data());

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

}

/*
 * Teil-C
 * Generate a triangle strip with N segments (2*N triangles)
 * in the xy plane (z=0).
 *
 * The vertices of the triangles should be in counter clock-wise order.
 *
 * The triangle strip must fill exactly the square [0, 1]x[0, 1].
 *
 * The first vertex should be at position (0,1,0) and the last
 * vertex at position (1,0,0)
 *
 * An example for N = 3:
 *
 *   ^
 *   |  ----------
 *   |  | /| /| /|
 * y |  |/ |/ |/ |
 *   |  ----------
 *   |
 *   |-------------->
 *           x
 *
 */
void generate_strip(
	std::uint32_t N,
	std::vector<glm::vec3>* vertices)
{
	cg_assert(N >= 1);
	cg_assert(vertices);

	vertices->clear();

	std::uint32_t num_vertices = 2 * (N + 1);
	float stride = 1.0 / N;

	for (int i = 0; i < num_vertices; i++)
	{
		if (i % 2 == 0)
		{
			vertices->push_back(glm::vec3(stride * i / 2.0,1, 0));
		}
		else
		{
			vertices->push_back(glm::vec3(stride * (i-1) / 2.0, 0, 0));
		}
	}

}

/*
 * Draw the given vertices as a triangle strip.
 * For each vertex, also set the corresponding color.
 */
void draw_triangle_strip(
	std::vector<glm::vec3> const& vertices,
	std::vector<glm::vec3> const& colors)
{
	cg_assert(vertices.size() == colors.size());

	int num_vertices = vertices.size();
	int n = int(num_vertices / 2 - 1);

	glBegin(GL_TRIANGLE_STRIP);
	for (int i = 0; i < vertices.size(); i++)
	{
		glColor3f(colors[i][0], colors[i][1], colors[i][2]);
		glVertex3f(vertices[i][0], vertices[i][1], vertices[i][2]);
	}
	glEnd();

}

/*
 * Integrate the given piecewise linear function
 * using trapezoidal integration.
 * The function is given at points
 *     x[0], ..., x[N]
 * and its corresponding values are
 *     y[0], ..., y[N]
 */
float integrate_trapezoidal(
	std::vector<float> const& x,
	std::vector<float> const& y)
{
	cg_assert(x.size() == y.size());
	cg_assert(x.size() > 1);

	int index = x.size();
	float T = 0.0;

	for (int i = 0; i < index - 1; i++)
	{
		T += 0.5 * (y[i] + y[i + 1])*(x[i + 1] - x[i]);
	}

	return T;
}

/*
 * Convert the given spectrum to RGB using your
 * implementation of integrate_trapezoidal(...)
 * The color matching functions and the wavelengths
 * for which they are given can be found in
 *     cglib/colors/cmf.h
 * and
 *     cglib/src/colors/cmf.cpp
 *
 * The wavelengths corresponding to the spectral values
 * given in spectrum are defined in cmf::wavelengths
 */
glm::vec3 spectrum_to_rgb(std::vector<float> const& spectrum)
{
	cg_assert(spectrum.size() == cmf::wavelengths.size());

	std::vector<float> x_s_array;
	std::vector<float> y_s_array;
	std::vector<float> z_s_array;

	for (int k = 0; k < spectrum.size(); k++)
	{
		x_s_array.push_back(spectrum[k] * cmf::x[k]);
		y_s_array.push_back(spectrum[k] * cmf::y[k]);
		z_s_array.push_back(spectrum[k] * cmf::z[k]);
	}

	float X = integrate_trapezoidal(cmf::wavelengths, x_s_array);
	float Y = integrate_trapezoidal(cmf::wavelengths, y_s_array);
	float Z = integrate_trapezoidal(cmf::wavelengths, z_s_array);

	glm::vec3 to_rgb = convert::xyz_to_rgb(glm::vec3(X, Y, Z));

	return to_rgb;
}
// CG_REVISION 719c9337e18b5a8d331f5b5580d53f2438ab503f
