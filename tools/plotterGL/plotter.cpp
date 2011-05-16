#include <iostream>
#include <vector>

#include <GL/glew.h> 
#include <GL/glut.h>

#include "plotter.h"

using namespace std;


const float cylHeight = 2.0f;
const float cylRadius = 5.0f;

float Plotter::angle[2] = {0.0f, 0.0f};
float Plotter::translate[3] = {0.0f, 0.0f, 0.0f};
Plotter* Plotter::instance = 0;
int Plotter::size = 0;
int Plotter::mouse_x = 0;
int Plotter::mouse_y = 0;
int Plotter::mouse_btn = -1;
bool Plotter::paintDisks = false;
bool Plotter::paintCylinder = true;

Plotter::Plotter() {
	if (instance != 0)
		throw "Only one Plotter instance a time allowed!";

	instance = this;
}

void Plotter::changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;

	float ratio =  w * 1.0 / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45.0f, ratio, 0.1f, 100.0f);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

int Plotter::loadData() {
	clog << "Loading data..." << endl;
	
	// prepare OpenGL vertex buffer
	GLuint vertexBufferID;
	glGenBuffers(1, &vertexBufferID);
	glBindBuffer(GL_ARRAY_BUFFER, vertexBufferID);

	// Load data into vector
	vector<float> data; // vector for data
	while (!cin.eof()) {
		float r;
		cin >> r; // Read over time value
		for (int i = 0; i < 3; i++) {
			cin >> r;
			data.push_back(r);
		}
		cin.ignore(1024, '\n'); // Skip over rest of line
	}

	// copy into array
	GLfloat *vertices = new GLfloat[data.size()];
	for (unsigned int i = 0; i < data.size(); i++) {
		vertices[i] = data[i];
	}
	clog << endl << "Got " << data.size() << " values." << endl;

	// copy into buffer
	glBufferData(GL_ARRAY_BUFFER, data.size()*sizeof(GLfloat),
			vertices, GL_STATIC_DRAW);
	clog << "Copied data into buffer." << endl;

	return data.size() - 1;
}

void Plotter::renderScene(void) {
	// Clear Color and Depth Buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();

	// Set the camera
	gluLookAt(	0.0f, 0.0f, 20.0f,
				0.0f, 0.0f,  0.0f,
				0.0f, 1.0f,  0.0f);

	glTranslatef(translate[0], translate[1], translate[2]);
	glRotatef(angle[0], 1.0f, 0.0f, 0.0f);
	glRotatef(angle[1], 0.0f, 1.0f, 0.0f);
	glTranslatef(0.0f, 0.0f, -cylHeight / 2.0f); // half height of model

	// set vertex pointer to 0 in buffer
	glVertexPointer(3, GL_FLOAT, 0, (GLvoid*)((char*) NULL));
	glColor4ub(255, 255, 255, 255);
	glDrawArrays(GL_LINE_STRIP, 0, size / 3);

	GLUquadric *quadric = gluNewQuadric();
	if (paintCylinder) {
		glColor4ub(255, 0, 0, 128);
		gluCylinder(quadric, cylRadius, cylRadius, cylHeight, 50, 50);
	}

	if (paintDisks) {
		glColor4ub(0, 255, 0, 128);
		gluDisk(quadric, 0, cylRadius, 50, 50);
		glTranslatef(0, 0, cylHeight);
		gluDisk(quadric, 0, cylRadius, 50, 50);
	}

	glutSwapBuffers();
}


void Plotter::init(int argc, char **argv, int w, int h) {
	// init GLUT and create window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	//glutInitWindowPosition(100,100);
	glutInitWindowSize(w, h);
	glutCreateWindow("PlotterGL");
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glEnable(GL_POLYGON_STIPPLE);

	// init GLEW
	glewInit();

	// load data from stdin
	size = loadData();

	// enable use of vertex arrays
	glEnableClientState(GL_VERTEX_ARRAY);

	// register callbacks
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	// glutIdleFunc(renderScene);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMotion);
	glutSpecialFunc(specialKey);

	// enter GLUT event processing cycle
	glutMainLoop();
}

void Plotter::specialKey(int key, int x, int y) {
	const float d = 0.1;

	switch (key) {
		case GLUT_KEY_UP:
			translate[1] += d;
			break;
		case GLUT_KEY_DOWN:
			translate[1] -= d;
			break;
		case GLUT_KEY_RIGHT:
			translate[0] += d;
			break;
		case GLUT_KEY_LEFT:
			translate[0] -= d;
			break;
		case GLUT_KEY_F2:
			paintDisks = !paintDisks;
			break;
		case GLUT_KEY_F1:
			paintCylinder = !paintCylinder;
			break;
		case GLUT_KEY_HOME:
			for (int i = 0; i < 3; i++)
				translate[i] = 0;
			for (int i = 0; i < 2; i++)
				angle[i] = 0;
			break;
	}

	glutPostRedisplay();
}

void Plotter::mouseButton(int btn, int state, int x, int y) {
	if (state == GLUT_UP) {
		mouse_btn = -1;
	}
	else {
		mouse_x = x;
		mouse_y = y;
		mouse_btn = btn;
	}

	glutPostRedisplay();
}

void Plotter::mouseMotion(int x, int y) {
	int dx = x - mouse_x;
	int dy = y - mouse_y;

	mouse_x = x;
	mouse_y = y;

	const float rotateSpeed = 0.5;
	const float translateSpeed = 0.1;

	if (mouse_btn == GLUT_LEFT_BUTTON) {
		angle[0] += (dy * rotateSpeed);
		angle[1] += (dx * rotateSpeed);
	}
	else if (mouse_btn == GLUT_RIGHT_BUTTON) {
		translate[2] += dy * translateSpeed;
	}

	glutPostRedisplay();
}

int main(int argc, char **argv) {
	Plotter p;
	p.init(argc, argv, 1000, 1000);
}
