/*
	Notes: Hermite Interpolation used for spline path of disc animation
	       Pressing "s" solves the puzzle from the initial state, not from any state
		   Rotation is done by changing orientation i.e 
										normal vector n = normalized(new_interpolated_position - previous_position)
										then finding out theta and rotating the disc with the angle
		   Disks do not penetrate through each other or the rods at any time.
		   Runs on Windows, Compile with Visual Studio. For Linux, change glut to freeglut when installing & including
*/

#include <iostream>
#include <glut.h>
#include <cmath>
#include <list>

using namespace std;

struct Vector3 {
	double x, y, z;
	Vector3() { x = y = z = 0.0; }
	Vector3(double x, double y, double z) : x(x), y(y), z(z) { }
	Vector3(Vector3 const& rhs) { *this = rhs; }
	Vector3& operator= (Vector3 const& rhs)
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		return *this;
	}
};

struct Disc {
	Disc() { normal = Vector3(0.0, 0.0, 1.0); }
	
	Vector3 position; //location
	Vector3 normal;   //orientation
};

struct ActiveDisc {    //Active Disc to be moved [later in motion]
	int disc_index;
	Vector3 start_pos, dest_pos;
	double u;		    // u E [0, 1]
	double step_u;
	bool is_in_motion;
	int direction;     // +1 for Left to Right & -1 for Right to left, 0 = stationary
};


// Rods and Discs Globals - Can be changed for different levels
const size_t NUM_DISCS = 4;
const double ROD_HEIGHT = 2.0;

struct Rod {
	Vector3 positions[NUM_DISCS];
	int occupancy_val[NUM_DISCS];
};

struct GameBoard {
	double x_min, y_min, x_max, y_max; //Base in XY-Plane
	double rod_base_rad;               //Rod's base radius
	Rod rods[3];	
};

struct solution_pair {
	size_t f, t;         //f = from, t = to
};

//Game Globals
Disc discs[NUM_DISCS];
GameBoard t_board;
ActiveDisc active_disc;
list<solution_pair> sol;
bool to_solve = false;

//Globals for window, time, FPS, moves
const size_t FPS = 30;
size_t moves = 0;
size_t prev_time = 0; 
size_t window_width = 600, window_height = 600;

const double PI = 3.141592653589793;

void initialize();
void initialize_game();
void display_handler();
void reshape_handler(int w, int h);
void keyboard_handler(unsigned char key, int x, int y);
void mouse_handler(int button, int state, int x, int y);
void anim_handler();
void move_disc(int from_rod, int to_rod);
Vector3 get_inerpolated_coordinate(Vector3 v1, Vector3 v2, double u);
void move_stack(int n, int f, int t);

int main(int argc, char** argv)
{

    // Initialize GLUT Window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(window_width, window_height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Towers of Hanoi");

	initialize();
	cout << "Towers of Hanoi" << endl;
	cout << "Press H for Help" << endl;

	//Callbacks
	glutDisplayFunc(display_handler);
	glutReshapeFunc(reshape_handler);
	glutKeyboardFunc(keyboard_handler);
	glutMouseFunc(mouse_handler);
	glutIdleFunc(anim_handler);

	glutMainLoop();
	return 0;
}

void initialize()
{
	glClearColor(0.0, 0.0, 0.0, 0.0); //Setting the clear color
	glShadeModel(GL_SMOOTH);		  //SMOOTH Shading

	glEnable(GL_DEPTH_TEST);		  //Enabling Depth Test

	//Setting Light0 parameters
	GLfloat light0_pos[] = { 0.0f, 0.0f, 0.0f, 1.0f }; // A positional light
	glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);

	glEnable(GL_LIGHTING);			  //Enabling Lighting
	glEnable(GL_LIGHT0);		      //Enabling Light0	


	//Globals initializations
	prev_time = glutGet(GLUT_ELAPSED_TIME);
	
	//Initializing Game State
	initialize_game();

}

void initialize_game()
{
	//Initializing 1)GameBoard t_board 2) Discs discs  3) ActiveDisc active_disc
	     // State
	
	//1) Initializing GameBoard
	t_board.rod_base_rad = 1.0;
	t_board.x_min = 0.0;
	t_board.x_max = 6.6 * t_board.rod_base_rad;
	t_board.y_min = 0.0;
	t_board.y_max = 2.2 * t_board.rod_base_rad;

	double x_center = (t_board.x_max - t_board.x_min) / 2.0;
	double y_center = (t_board.y_max - t_board.y_min) / 2.0;

	double dx = (t_board.x_max - t_board.x_min) / 3.0; //Since 3 rods
	double r = t_board.rod_base_rad;

    //Initializing Rods Occupancy value		
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t h = 0; h < NUM_DISCS; h++)
		{
			if (i == 0) 
			{
				t_board.rods[i].occupancy_val[h] = NUM_DISCS - 1 - h;
			}
			else t_board.rods[i].occupancy_val[h] = -1;
		}
	}

	
	//Initializing Rod positions
	for (size_t i = 0; i < 3; i ++)
	{
		for (size_t h = 0; h < NUM_DISCS; h++)
		{
			double x = x_center + ((int)i - 1) * dx;
			double y = y_center;
			double z = (h + 1) * 0.3;
			Vector3& pos_to_set = t_board.rods[i].positions[h];
			pos_to_set.x = x;
			pos_to_set.y = y;
			pos_to_set.z = z;	
		}
	}


	//2) Initializing Discs
	for (size_t i = 0; i < NUM_DISCS; i++)
	{
		discs[i].position = t_board.rods[0].positions[NUM_DISCS - i - 1];
		//Normals are initialized whie creating a Disc object - ie in constructor of Disc
	}


	//3) Initializing Active Disc
	active_disc.disc_index = -1;
	active_disc.is_in_motion = false;
	active_disc.step_u = 0.02;         //Verify
	active_disc.u = 0.0;
	active_disc.direction = 0;
}

//Draw function for drawing a cylinder given position and radius and height
void draw_solid_cylinder(double x, double y, double r, double h)
{
	GLUquadric* q = gluNewQuadric();
	GLint slices = 50;
	GLint stacks = 10;
	
	glPushMatrix();
	glTranslatef(x, y, 0.0f);
	gluCylinder(q, r, r, h, slices, stacks);
	glTranslatef(0, 0, h);
	gluDisk(q, 0, r, slices, stacks);
	glPopMatrix();
	gluDeleteQuadric(q);
}

//Draw function for drawing rods on a given game board i.e. base
void draw_board_and_rods(GameBoard const& board)
{
	//Materials, 
	GLfloat mat_white[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat mat_yellow[] = { 1.0f, 1.0f, 0.0f, 1.0f };

	glPushMatrix();
	//Drawing the Base Rectangle [where the rods are placed]
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_white);
	glBegin(GL_QUADS);
		glNormal3f(0.0f, 0.0f, 1.0f);
		glVertex2f(board.x_min, board.y_min);
		glVertex2f(board.x_min, board.y_max);
		glVertex2f(board.x_max, board.y_max);
		glVertex2f(board.x_max, board.y_min);
	glEnd();

	//Drawing Rods and Pedestals
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_yellow);

	double r = board.rod_base_rad;
	for (size_t i = 0; i < 3; i++)
	{
		Vector3 const& p = board.rods[i].positions[0];
		draw_solid_cylinder(p.x, p.y, r * 0.1, ROD_HEIGHT - 0.1);
		draw_solid_cylinder(p.x, p.y, r, 0.1);
	}

	glPopMatrix();
}


// Draw function for drawing discs
void draw_discs()
{
	int slices = 100;
	int stacks = 10;

	double rad;

	GLfloat r, g, b;
	GLfloat emission[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat no_emission[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat material[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	for (size_t i = 0; i < NUM_DISCS; i++)
	{
		switch (i)
		{
		case 0: r = 1.0f; g = 0.0f; b = 0.0f;
			break;
		case 1: r = 0.0f; g = 1.0f; b = 0.0f;
			break;
		case 2: r = 0.0f, g = 0.0f; b = 1.0f;
			break;
		default: r = g = b = 1.0f;
			break;
		};

		material[0] = r;
		material[1] = g;
		material[2] = b;
		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material);

		GLfloat u = 0.0f;
		
		if (i == active_disc.disc_index) {
			glMaterialfv(GL_FRONT, GL_EMISSION, emission);
			u = active_disc.u;			
		}

		GLfloat factor = 1.0f;
		switch (i) {
			case 0: factor = 0.2;
				break;
			case 1: factor = 0.4;
				break;
			case 2: factor = 0.6;
				break;
			case 3: factor = 0.8;
				break;
			default: break;
		};
		rad = factor * t_board.rod_base_rad;
		int d = active_disc.direction;


		glPushMatrix();
		glTranslatef(discs[i].position.x, discs[i].position.y, discs[i].position.z);
		double theta = acos(discs[i].normal.z);
		theta *= 180.0f / PI;
		glRotatef(d * theta , 0.0f, 1.0f, 0.0f);
		glutSolidTorus(0.2 * t_board.rod_base_rad, rad, stacks, slices);
		glPopMatrix();

		glMaterialfv(GL_FRONT, GL_EMISSION, no_emission); 
	}
}

void display_handler()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	double x_center = (t_board.x_max - t_board.x_min) / 2.0;
	double y_center = (t_board.y_max - t_board.y_min) / 2.0;
	double r = t_board.rod_base_rad;

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(	x_center, y_center - 10.0, 3.0 * r,
				x_center, y_center, 3.0,
				0.0, 0.0, 1.0  );

	glPushMatrix();
		draw_board_and_rods(t_board);
		draw_discs();	
	glPopMatrix();
	glFlush();
	glutSwapBuffers();
}

void reshape_handler(int w, int h)
{
	window_width = w;
	window_height = h;

	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (GLfloat)w / (GLfloat)h, 0.1, 20.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void move_stack(int n, int f, int t)
{
	if (n == 1) {
		solution_pair s;
		s.f = f;
		s.t = t;
		sol.push_back(s);         //pushing the (from, to) pair of solution to a list [so that it can be animated later]
		moves++;
		cout << "From rod " << f << " to Rod " << t << endl;
		return;
	}
	move_stack(n - 1, f, 3 - t - f);
	move_stack(1, f, t);
	move_stack(n - 1, 3 - t - f, t);
}

//Solve from 1st rod to 2nd
void solve()
{
	move_stack(NUM_DISCS, 0, 2);
}

void keyboard_handler(unsigned char key, int x, int y)
{

    //Console Outputs
	switch (key) 
	{
	case 27: cout << "Steps Taken: " << moves << endl;
		exit(0);
		break;

	case 'h':
	case 'H': cout << "ESC: Quit" << endl;
		cout << "S: Solve from Initial State" << endl;
		cout << "Left Click: to select/deselect the disc and to specify destination rod" << endl;
		cout << "T: shows the current number of steps used" << endl;
		cout << "H: Help" << endl;
		break;

	case 's':
	case 'S': 
		if (t_board.rods[0].occupancy_val[NUM_DISCS - 1] < 0)
			break;
		
		solve();
		to_solve = true;
		break;

	case 't':
	case 'T':
		cout << "Steps = " << moves << endl;
		break;

	default:
		break;
	};
}

void mouse_handler(int button, int state, int x, int y)
{
	if (active_disc.is_in_motion) return;

	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		size_t rod = (int) x / (window_width / 3.0);
		//std::cout << rod << std::endl;

		int i;
		for (i = NUM_DISCS - 1; i >=0 && t_board.rods[rod].occupancy_val[i] < 0; i--);

		int selected_disc_i = -1;
		//cout << selected_disc_i << endl;
		if ( (i >= 0) && (t_board.rods[rod].occupancy_val[i] >= 0) )
		{
			selected_disc_i = t_board.rods[rod].occupancy_val[i];
		}
				
		//cout << selected_disc_i << endl;

		if (active_disc.disc_index < 0 && selected_disc_i >=0 ) {
			active_disc.disc_index = selected_disc_i;
		}

		else if (active_disc.disc_index == selected_disc_i) {
			active_disc.disc_index = -1;
		}

		else {
			//cout << "Yes" << endl;
			if (active_disc.disc_index < 0) return;
			if (selected_disc_i >= 0 && active_disc.disc_index > selected_disc_i) return;


			size_t to_rod = rod, from_rod;
			for (size_t i = 0; i < 3; i++) {
				for (size_t j = 0; j < NUM_DISCS; j++) {
					if (t_board.rods[i].occupancy_val[j] == active_disc.disc_index)
						from_rod = i;
				}
			}

			move_disc(from_rod, to_rod);
			moves++; //Steps++
		}
	}
	glutPostRedisplay();
}

void move_disc(int from_rod, int to_rod)
{

	int d = to_rod - from_rod;
	if (d > 0) active_disc.direction = 1;
	else if (d < 0) active_disc.direction = -1;

	if ((from_rod == to_rod	) || (from_rod < 0) || (to_rod < 0) || (from_rod > 2) || (to_rod > 2) )
		return;
		
	size_t i;
	for (i = NUM_DISCS - 1; i >= 0 && t_board.rods[from_rod].occupancy_val[i] < 0; i--);
	if ( (i < 0) || (i == 0 && t_board.rods[from_rod].occupancy_val[i] < 0) )
		return; //Either the index < 0 or index at 0 and occupancy < 0 => it's an empty rod

	active_disc.start_pos = t_board.rods[from_rod].positions[i];
	
	active_disc.disc_index = t_board.rods[from_rod].occupancy_val[i];
	active_disc.is_in_motion = true;
	active_disc.u = 0.0;


	size_t j;
	for (j = 0; j < NUM_DISCS - 1 && t_board.rods[to_rod].occupancy_val[j] >= 0; j++);
	active_disc.dest_pos = t_board.rods[to_rod].positions[j];
	
	t_board.rods[from_rod].occupancy_val[i] = -1;
	t_board.rods[to_rod].occupancy_val[j] = active_disc.disc_index;
}


Vector3 get_inerpolated_coordinate(Vector3 sp, Vector3 tp, double u)
{
	//4 Control points
	Vector3 p;
	double x_center = (t_board.x_max - t_board.x_min) / 2.0;
	double y_center = (t_board.y_max - t_board.y_min) / 2.0;

	double u3 = u * u * u;
	double u2 = u * u;

	Vector3 cps[4]; //P1, P2, dP1, dP2

	
	//Hermite Interpolation [Check Reference for equation of spline]
	{
		//P1
		cps[0].x = sp.x;
		cps[0].y = y_center;
		cps[0].z = ROD_HEIGHT + 0.2 * (t_board.rod_base_rad);

		//P2
		cps[1].x = tp.x;
		cps[1].y = y_center;
		cps[1].z = ROD_HEIGHT + 0.2 * (t_board.rod_base_rad);

		//dP1
		cps[2].x = (sp.x + tp.x) / 2.0 - sp.x ;
		cps[2].y = y_center;
		cps[2].z = 2 * cps[1].z; //change 2 * ..

		//dP2
		cps[3].x = tp.x - (tp.x + sp.x)/2.0;
		cps[3].y = y_center;
		cps[3].z = -cps[2].z; //- cps[2].z;


		double h0 = 2 * u3 - 3 * u2 + 1;
		double h1 = -2 *u3 + 3 * u2;
		double h2 = u3 - 2 * u2 + u;
		double h3 = u3 - u2;

		p.x = h0 * cps[0].x + h1 * cps[1].x + h2 * cps[2].x +  h3 * cps[3].x;
		p.y = h0 * cps[0].y + h1 * cps[1].y + h2 * cps[2].y +  h3 * cps[3].y;
		p.z = h0 * cps[0].z + h1 * cps[1].z + h2 * cps[2].z +  h3 * cps[3].z;
	
	}

	return p;
}

//Normalize function for a vector
void normalize(Vector3& v)
{
	double length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	if (length == 0.0) return;
	v.x /= length;
	v.y /= length;
	v.z /= length;
}


Vector3 operator-(Vector3 const& v1, Vector3 const& v2)
{
	return Vector3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}



void anim_handler()
{
	int curr_time = glutGet(GLUT_ELAPSED_TIME);
	int elapsed = curr_time - prev_time; // in ms
	if (elapsed < 1000 / FPS) return;

	prev_time = curr_time;

	if (to_solve && active_disc.is_in_motion == false) {
		solution_pair s = sol.front();
		
        cout << s.f << ", " << s.t << endl;
		
		sol.pop_front();
		int i;
		for (i = NUM_DISCS; i >= 0 && t_board.rods[s.f].occupancy_val[i] < 0; i--);
		int ind = t_board.rods[s.f].occupancy_val[i];
		
		if (ind >= 0)
		active_disc.disc_index = ind;
		
		move_disc(s.f, s.t);
		if (sol.size() == 0)
			to_solve = false;
	}

	if (active_disc.is_in_motion)
	{
		int ind = active_disc.disc_index;
		ActiveDisc& ad = active_disc;

		if (ad.u == 0.0 && (discs[ind].position.z < ROD_HEIGHT + 0.2 * (t_board.rod_base_rad)) )
		{
			discs[ind].position.z += 0.05;
			glutPostRedisplay();
			return;
		}

		static bool done = false;
		if (ad.u == 1.0 && discs[ind].position.z > ad.dest_pos.z )
		{
			done = true;
			discs[ind].normal = Vector3(0, 0, 1);
			discs[ind].position.z -= 0.05;
			glutPostRedisplay();
			return;
		}

		ad.u += ad.step_u;
		if (ad.u > 1.0) {
			ad.u = 1.0;			
		}

		if (!done) {
			Vector3 prev_p = discs[ind].position;
			Vector3 p = get_inerpolated_coordinate(ad.start_pos, ad.dest_pos, ad.u);
			discs[ind].position = p;
			discs[ind].normal.x = (p - prev_p).x;
			discs[ind].normal.y = (p - prev_p).y;
			discs[ind].normal.z = (p - prev_p).z;
			normalize(discs[ind].normal);
		}

		if (ad.u >= 1.0 && discs[ind].position.z <= ad.dest_pos.z) {
			discs[ind].position.z = ad.dest_pos.z;
			ad.is_in_motion = false;
			done = false;
			ad.u = 0.0;
			discs[ad.disc_index].normal = Vector3(0, 0, 1);
			ad.disc_index = -1;
			
		}
		glutPostRedisplay();
	}
}