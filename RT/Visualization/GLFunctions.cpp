#ifdef _WIN32 || _WIN64
	//#pragma comment(lib, "pthreadVSE1.lib")
    #include <windows.h>
    #include <GL/gl.h>
    #include <GL/freeglut.h>
#include <GL/glu.h>
#elif __linux__   
    #include <GL/gl.h>
    #include <GL/glut.h>
    #include <GL/glu.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath>

#include <time.h>

#include "GLFunctions.h"
#include "math_and_physical_constants.h"
#include "ray_path.h"

//#include <pthread.h>

int keyb[256]; //keyboard list

std::vector<vec3> mesh_colors;

std::vector<vec3> ray_colors;

std::vector<std::vector<inter_point>> paths;		// --- Each element of this array is a path. Each path is a list of points whose coordinates are stored
													//     in the post field.

extern std::vector<path> global_path_list;

extern std::vector<TX> txs;
extern std::vector<RX> rxs;

/***********************/
/* INPUT LOOP FUNCTION */
/***********************/
void input_loop() {

	paths.clear();

	for (size_t i = 0; i < global_path_list.size(); i++) {
   
		std::vector<inter_point> path;
   
		for(size_t k = 0; k < global_path_list[i].intersections.size(); k++) {

			inter_point p = { global_path_list[i].intersections[k].type, global_path_list[i].intersections[k].pos };
			path.push_back(p);

		}

		paths.push_back(path);
	
	}

}

/*****************/
/* CAMERA STRUCT */
/*****************/
struct camera{
  vec3    pos;
  GLfloat elevation;
  GLfloat azimut;
  GLfloat zangle;
}cam;

/*********************************/
/* DRAW REFERENCE FRAME FUNCTION */
/*********************************/
void draw_reference_frame(){

	glLineWidth(1.5f);
	glBegin(GL_LINES);

	// --- x-axis
	glColor3f (1.0f,  .0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(1.0f, 0.0f, 0.0f);

	// --- y-axis
	glColor3f (0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 1.0f, 0.0f);

	// --- z-axis
	glColor3f (0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 1.0f);

	glEnd();
	glLineWidth(1.0f);
}

/***********************/
/* KEYBOARDUP FUNCTION */
/***********************/
void keyboardUp(unsigned char key, int x, int y){

    //bool show_index = false;
    keyb[key] = 0;

//     if(key == 'v'){
//         if(!current->is_leaf()){
//             current = reinterpret_cast<inner_node*>(current)->nodes_ptr[0];
//             show_index = true;
//         }
//     }
//     if(key ==  'b'){
//         if(!current->is_leaf()){
//             current = reinterpret_cast<inner_node*>(current)->nodes_ptr[1];
//             show_index = true;
//         }
//     }
//     if(key ==  'p'){
//         if(current->parent != NULL){
//             current = current->parent;
//             show_index = true;
//         }
//     }
//
//     if( show_index && !current->is_leaf() ){
//         std::cout << "node1 idx:" <<  reinterpret_cast<inner_node*>(current)->nodes_idx[0]
//                   << "node2 idx:" <<  reinterpret_cast<inner_node*>(current)->nodes_idx[1]
//                   << std::endl;
//     }

}

/*************************/
/* KEYBOARDDOWN FUNCTION */
/*************************/
void keyboardDown(unsigned char key, int x, int y){

  keyb[key] = 1;
}

/*****************/
/* LOOP FUNCTION */
/*****************/
void loop(int value){
    GLfloat dx,dz;
    if(keyb['w']){
        dx = sinf(cam.azimut*M_PI/180.0f)/10.0f;
        dz = -cosf(cam.azimut*M_PI/180.0f)/10.0f;
        cam.pos.x +=dx;
        cam.pos.z +=dz;
    }
    if(keyb['s']){
        dx = sinf(cam.azimut*M_PI/180.0f)/10.0f;
        dz = -cosf(cam.azimut*M_PI/180.0f)/10.0f;
        cam.pos.x -= dx;
        cam.pos.z -= dz;
    }

    if(keyb['a']){
        dx = -cosf(cam.azimut*M_PI/180.0f)/10.0f;
        dz = -sinf(cam.azimut*M_PI/180.0f)/10.0f;
        cam.pos.x +=dx;
        cam.pos.z +=dz;
    }
    if(keyb['d']){
        dx = -cosf(cam.azimut*M_PI/180.0f)/10.0f;
        dz = -sinf(cam.azimut*M_PI/180.0f)/10.0f;

        cam.pos.x -= dx;
        cam.pos.z -= dz;
    }


    if(keyb['m']){
        cam.pos.y+=0.1;
    }
    if(keyb['n']){
        cam.pos.y-=0.1;
    }
    if(keyb['i']){
        cam.elevation+=0.5f;
    }
    if(keyb['k']){
        cam.elevation-=0.5f;
    }
    if(keyb['j']){
        cam.azimut-=0.5f;
    }
    if(keyb['l']){
        cam.azimut+=0.5f;
    }


    if(keyb[27]==1){
        exit(0);
    }
    glutPostRedisplay();
    glutTimerFunc(30,loop, 0);
}

/*****************/
/* INIT FUNCTION */
/*****************/
void Visualize(int *argc, char** argv) {
	
	glutInit(argc, argv);											// --- Initializes the GLUT library
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);		// --- Sets the initial display mode
    //glutInitWindowSize(1024, 768);									// --- Sets the initial window size
    glutInitWindowSize(2048, 1024);									// --- Sets the initial window size
    glutInitWindowPosition (100, 100);								// --- Sets the initial window position
    glutCreateWindow(argv[0]);										// --- Creates a top level window

	// --- Black background
	//glClearColor(0.0f, 0.0f, 0.0f, 0.0f);							// --- Clears the color buffer
	// --- White background
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);							// --- Clears the color buffer
    glClearDepth(1.0f);												// --- Specifies the depth value used when the depth buffer is cleared
    glDepthFunc(GL_LESS);											// --- Specifies the value used for depth buffer comparisons
    glEnable(GL_DEPTH_TEST);										// --- Specifies a symbolic constant indicating a GL capability
    glShadeModel(GL_SMOOTH);										// --- Specifies a symbolic value representing a shading technique

	srand((unsigned int)time(NULL));

	cam.pos = make_vec3(0, 3, 15);
    cam.elevation = 0;
    cam.azimut = 0;

	input_loop();
	//pthread_t thread1;

 //   int  iret1 = pthread_create(&thread1, NULL, input_loop, NULL);
 //   if(iret1) {
 //       fprintf(stderr,"Error - pthread_create() return code: %d\n",iret1);
 //       exit(EXIT_FAILURE);
 //   }

    mesh_colors.push_back(make_vec3(0.5,0,0));
    mesh_colors.push_back(make_vec3(0,0.5,0));
    mesh_colors.push_back(make_vec3(0,0,0.5));
    mesh_colors.push_back(make_vec3(1,1,0));
    mesh_colors.push_back(make_vec3(0,1,1));
    mesh_colors.push_back(make_vec3(1,0,1));
    mesh_colors.push_back(make_vec3(1,0.5,0.5));
    mesh_colors.push_back(make_vec3(0.2,0.5,1));
    mesh_colors.push_back(make_vec3(0.2,0.2,0.2));

    ray_colors.push_back(make_vec3(0.2f, 0.8f, 1.0f));
    ray_colors.push_back(make_vec3(0.2f, 1.0f, 1.0f));
    ray_colors.push_back(make_vec3(0.2f, 1.0f, 1.0f));
    ray_colors.push_back(make_vec3(1.0f, 0.2f, 0.2f));
    ray_colors.push_back(make_vec3(0.0f, 0.2f, 1.0f));

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboardDown);
    glutKeyboardUpFunc(keyboardUp);
    glutTimerFunc(50, loop, 0);
    glutMainLoop();
}

/**************************/
/* DRAW RECEIVER FUNCTION */
/**************************/
void draw_rx(const RX &rx){

	glPushMatrix();
    glLineWidth(1.5f);
    glBegin(GL_LINES);
    glColor3f (1.0f, .0f, 0.0f);
    draw_segment(rx.pos,rx.pos+0.15*rx.x_v);

    glColor3f (0.0f, 1.0f, 0.0f);
    draw_segment(rx.pos,rx.pos+0.15*rx.y_v);

    glColor3f (0.0f, 0.0f, 1.0f);
    draw_segment(rx.pos-0.10*rx.z_v,rx.pos+0.20*rx.z_v);

    glEnd();
    glLineWidth(1.0f);
    glColor3f (1.0f, 1.0f, 0.0f);

    glTranslatef(rx.pos.x,rx.pos.y,rx.pos.z);
    glutSolidSphere(0.1, 20, 20);
    glPopMatrix();
}

/*****************************/
/* DRAW TRANSMITTER FUNCTION */
/*****************************/
void draw_tx(const TX &tx){

	glPushMatrix();
    glLineWidth(1.5f);
    glBegin(GL_LINES);
    glColor3f (1.0f, .0f, 0.0f);

    draw_segment(tx.pos, tx.pos + 0.15 * tx.x_v);

    glColor3f (0.0f, 1.0f, 0.0f);
    draw_segment(tx.pos, tx.pos + 0.15 * tx.y_v);

    glColor3f (0.0f, 0.0f, 1.0f);
    draw_segment(tx.pos - 0.10 * tx.z_v, tx.pos + 0.20 * tx.z_v);

    glEnd();
    glLineWidth(1.0f);
    glColor3f (1.0f, 0.0f, 1.0f);

    glTranslatef(tx.pos.x, tx.pos.y, tx.pos.z);
    //glutSolidSphere(0.05, 20, 20);
    glutSolidSphere(0.1, 20, 20);

    glPopMatrix();
}

/********************/
/* DISPLAY FUNCTION */
/********************/
void display(void){

	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    glPushMatrix();

    // --- Positioning the cam on azimuthal installation
    glRotatef(cam.elevation, 1.0, 0.0, 0.0);
    glRotatef(cam.azimut, 0.0, 1.0, 0.0);
    glTranslatef(-cam.pos.x, -cam.pos.y, -cam.pos.z);
    glRotatef(-90, 1.0, 0.0, 0.0);

    // --- Axes
    glPushMatrix();
    glScalef(1.0, 1.0, 1.0);
    draw_reference_frame();
    glPopMatrix();

    // --- Draw plane
    if(!mesh_loaded) {
    	glColor3f (1.0, 1.0, 0.0);
        draw_plane(10,10);
    } else draw_model();

    // --- Draw tx antennas
    for(size_t i = 0; i < txs.size(); i++) draw_tx(txs[i]);

    // --- Draw rx antennas
    for(size_t i = 0; i < rxs.size(); i++) draw_rx(rxs[i]);

    draw_rays();

    glPopMatrix();

    glutSwapBuffers();
}

/***********************/
/* DRAW MODEL FUNCTION */
/***********************/
void draw_model(void){

	if(mesh.verts.size() != 0){

		GLfloat linelen = mesh.min_len * 0.7f;
//         // activate and specify pointer to vertex array
//         glEnableClientState(GL_VERTEX_ARRAY);
//         glVertexPointer(3, GL_FLOAT, 0, mesh.vertices);
//
//         // draw a cube
//         glDrawElements(GL_TRIANGLES, 3*mesh.flist.size(), GL_UNSIGNED_INT, mesh.indices);
//
//         // deactivate vertex arrays after drawing
//         glDisableClientState(GL_VERTEX_ARRAY);

		glBegin(GL_TRIANGLES);

		for (size_t i = 0; i < mesh.faces.size(); i++) {

			vec3 color = mesh_colors[mesh.faces[i].w % mesh_colors.size()];

			glColor3f (color.x, color.y, color.z);

			float4  a = mesh.verts[mesh.faces[i].x];
			float4  b = mesh.verts[mesh.faces[i].y];
			float4  c = mesh.verts[mesh.faces[i].z];
         ////GLfloat ka = 1/(rho1[faces[i].x]*rho2[faces[i].x]);
         ////GLfloat kb = 1/(rho1[faces[i].y]*rho2[faces[i].y]);
         ////GLfloat kc = 1/(rho1[faces[i].z]*rho2[faces[i].z]);
         ////vec3f ca = get_colour(ka,min_k,max_k);
         ////vec3f cb = get_colour(kb,min_k,max_k);
         ////vec3f cc = get_colour(kc,min_k,max_k);
         //////printf("%f %f %f \n",ca.x, ca.y, ca.z);

         ////glColor3f (ca.x, ca.y, ca.z);
			glVertex3f(a.x, a.y, a.z);

             ////glColor3f (cb.x, cb.y, cb.z);
			glVertex3f(b.x, b.y, b.z);

         ////glColor3f (cc.x, cc.y, cc.z);
			glVertex3f(c.x, c.y, c.z);

		}

		glEnd();

		// --- draw normals
		glLineWidth(1.5f);

		glBegin(GL_LINES);
		for(size_t i =0; i< mesh.verts.size(); i++) {
			float4  a = mesh.verts[i];
			float4  n = mesh.normals_and_k1[i];
			glColor3f (0.5f, 1.0f, 1.0f);
			glVertex3f(a.x, a.y, a.z);
			glVertex3f(a.x + linelen * n.x, a.y + linelen * n.y, a.z + linelen * n.z);
		}
		glEnd();

		// --- draw_X1
		glBegin(GL_LINES);
		for(size_t i = 0; i < mesh.verts.size(); i++) {
			float4  a = mesh.verts[i];
			float4  x1 = mesh.x1_and_k2[i];
			glColor3f (0.1f, 1.0f, 0.3f);
			glVertex3f(a.x, a.y, a.z);
			glVertex3f(a.x + linelen * x1.x, a.y + linelen * x1.y, a.z + linelen * x1.z);
		}
		glEnd();

		glLineWidth(1.0f);

    }
}

/**********************/
/* DRAW RAYS FUNCTION */
/**********************/
void draw_rays(){
    glPushMatrix();

    for(size_t i = 0; i < paths.size(); i++){
        glBegin(GL_LINE_STRIP);
        for(size_t j = 0; j < paths[i].size(); j++) {
        	vec3 color = ray_colors[paths[i][j].type];
            glColor4f (color.x, color.y, color.z,0.4f);
            glVertex3f(paths[i][j].pos.x,paths[i][j].pos.y,paths[i][j].pos.z);
        }
        glEnd();
    }
    glPopMatrix();

}

/********************/
/* RESHAPE FUNCTION */
/********************/
void reshape (int w, int h){
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(65.0f,(GLfloat) w/((GLfloat)h),0.01f,1000.0f);
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();
}

/*************************/
/* DRAW SEGMENT FUNCTION */
/*************************/
void draw_segment(const vec3 A, const vec3 B){
    glVertex3f(A.x,A.y,A.z);
    glVertex3f(B.x,B.y,B.z);
}

/***********************/
/* DRAW PLANE FUNCTION */
/***********************/
void draw_plane(int Nx,int Ny){

	int i, j;
    GLfloat x, y;

    for(i = 0; i < 2 * Nx; i++) {
    	for(j = 0; j < 2 * Ny; j++) {

    		x = i - Nx + 0.5f;
            y = j - Ny + 0.5f;

            if((i % 2 == 0 && j % 2 == 0) || (i % 2 == 1 && j % 2 == 1)) glColor3f (1.0f, 1.0f, 1.0f);
            else glColor3f (0.3f,0.3f,0.4f);

            glBegin(GL_TRIANGLES);
            glVertex3f(x-0.5f, y+0.5f, 0);
            glVertex3f(x+0.5f, y-0.5f, 0);
            glVertex3f(x-0.5f, y-0.5f, 0);
            glVertex3f(x-0.5f, y+0.5f, 0);
            glVertex3f(x+0.5f, y+0.5f, 0);
            glVertex3f(x+0.5f, y-0.5f, 0);
            glEnd();
        }
    }
}
