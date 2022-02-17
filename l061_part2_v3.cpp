//Trisha Rayan
//Lab 6 Part 2 Canny Edge Detection
#include <bits/stdc++.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>   
#include <cstdlib> 
#include <bits/stdc++.h> 
#include <math.h> 
#include <stdlib.h>   
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

const int MAX_PIXEL = 800;
const int MAX_CLUSTER = 4;
const int MAX_POINTS = 100;
int PPM_WIDTH;
int PPM_HEIGHT;
const int THRES = 120;
const int LOW_THRES = 30;
const int HIGH_THRES = 100;
int MAX_VOTE = 0;
const int CP_THRES = 40; 
//int PPM_SIZE;
//int MAX_RGB = 255;



class Point
{
	private:
		double X, Y;
	
	public:
		//defualt constructor 
		Point () {X=0; Y=0;}
        Point (double a, double b) {X=a; Y=b;}
		
		//setter function
		void setPoint(double a, double b) { X = a; Y = b;}	
        void setX(double a) { this->X = a;}
        void setY(double b) {this->Y = b;}

		//getter functions
		double getX(void) { return X;} 
		double getY(void) {return Y;}
};

class Color
{
	private:
		int r, g, b;
	
	public:
		//defualt constructor 
		Color () {r=0; g=0; b=0;}
        Color (int red, int green, int blue) {r=red; g=green; b=blue;}
		
		//setter function
		void setColor(int red, int green, int blue) {r=red; g=green; b=blue;}	

        //getter functions
		int getR(void) {return r;}
		int getG(void) {return g;}
        int getB(void) {return b;}
};

class Line
{
	private:
		int x1, y1, x2, y2;
	
	public:
		//defualt constructor 
		Line () {x1=0; y1=0; x2=0; y2=0;}
        Line (int ax, int ay, int bx, int by) {x1=ax; y1=ay; x2=bx; y2=by;}
		
		//setter function
		void setPoint1(int ax, int ay) { x1 = ax; y1 = ay;}	
        void setPoint2(int ax, int ay) { x2 = ax; y2 = ay;}	



		//getter functions
		int getP1X(void) { return x1;} 
		int getP1Y(void) {return y1;}
        int getP2X(void) { return x2;} 
		int getP2Y(void) {return y2;}
};

//Color col[MAX_CLUSTER] {{255,0,0}, {0,0,255}, {239,0,255}, {0,255,0}, {0,205,255}};
Color col[MAX_CLUSTER] {{0,0,0}, {85,85,85}, {170,170,170}, {255,255,255}}; //init cluster

typedef struct Pixel {
    int r;
    int g;
    int b;

    Pixel() {r = 0; g = 0; b = 0;    }

    Pixel(int r, int g, int b) {
    this->r = r;
    this->g = g;
    this->b = b;
    }
}pix;

static pix ppm[MAX_PIXEL][MAX_PIXEL];

bool sortbysecdesc(const pair<string,int> &a, const pair<string,int> &b){
    return a.second > b.second;
}


double randval(double min, double max) {
 return (double) (min + (max - min) * ((double) rand() / RAND_MAX)); 
}

int randval(int min, int max) {
 return (int) (min + (max - min) * ((int) rand() / RAND_MAX)); 
}

list<double> generatePoints(int total)
{
    list<double> pntList;
    for(int i=0; i < total ; ++i)
    {
        double pnt = randval(0.0, 1.0);
        //cout << std::fixed << std::setprecision(23) << pnt << endl;
        pntList.push_back(pnt); 
    }  
    return pntList;
}

string convertToString(char* a) 
{ 
    string s(a);
    return s;
}

double distance(Point p1, Point p2)  
{  
    return sqrt( (p1.getX() - p2.getX())*(p1.getX() - p2.getX()) +  
                (p1.getY() - p2.getY())*(p1.getY() - p2.getY())  
            );  
}

double distanceRGB(Color c1, Color c2)  
{  
     //sqrt((l1 - l2)^2 + (a1 - a2)^2 + (b1 - b2)^2)
    return sqrt( (c1.getR() - c2.getR())*(c1.getR() - c2.getR()) +  
                (c1.getG() - c2.getG())*(c1.getG() - c2.getG())  +
                (c1.getB() - c2.getB())*(c1.getB() - c2.getB())
            );  
}

std::vector<Point> vect_p; 
std::vector<Color> vect_ppm_p; 
std::vector<Color> vect_m;

void drawPoint(vector<Color> &col, int x, int y, int r, int g, int b)
{
    if (x < 0 || x >= MAX_PIXEL || y < 0 || y >= MAX_PIXEL)
        return;
    col[y*PPM_WIDTH+x] = Color(r,g,b);
}

void createCircle(vector<Color> &col, int midx, int midy, double rad, int r, int g, int b)
{
    int x, y, xmax, y2, y2_new, ty;
    xmax = (int)(rad * 0.70710678) + 1;
    y = (int)rad;
    y2 = y * y;
    ty = (2 * y) - 1;
    y2_new = y2;

    for (x = 0; x <= xmax; x++)
    {
        if ((y2 - y2_new) >= ty)
        {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }
        drawPoint(col, x + midx, y + midy, r, g, b);
        drawPoint(col,x + midx, -y + midy, r, g, b);
        drawPoint(col,-x + midx, y + midy, r, g, b);
        drawPoint(col,-x + midx, -y + midy, r, g, b);
        drawPoint(col,y + midx, x + midy, r, g, b);
        drawPoint(col,y + midx, -x + midy, r, g, b);
        drawPoint(col,-y + midx, x + midy, r, g, b);
        drawPoint(col,-y + midx, -x + midy, r, g, b);

        y2_new -= (2 * x) - 3;
    }
}

int **promoteWeakToStrong(int **hysteresis, int r, int c, bool **visited)
{
        //cout << endl << "find:" << r << ":" << c << ":" << arr[r][c];
    if((r<0) or (c<0) or (r>= PPM_HEIGHT) or (c>=PPM_WIDTH) or (hysteresis[r][c] == 0) or (visited[r][c] == true))
         return hysteresis;

    //cout << endl << "find:" << r << ":" << c << ":" << hysteresis[r][c];
    visited[r][c]=true;
    if(hysteresis[r][c] == 1)
    {
        hysteresis[r][c] = 2; 
        //cout << endl << "replaced:" << r << ":" << c << ":" << hysteresis[r][c];
    } 

    //arr[r][c] = 2;
    //int num = 0;
    promoteWeakToStrong(hysteresis,r, c+1, visited);
    promoteWeakToStrong(hysteresis,r-1, c, visited);
    promoteWeakToStrong(hysteresis,r, c-1, visited);
    promoteWeakToStrong(hysteresis,r+1, c, visited);
    promoteWeakToStrong(hysteresis,r-1, c+1, visited);
    promoteWeakToStrong(hysteresis,r-1, c-1, visited);
    promoteWeakToStrong(hysteresis,r+1, c-1, visited);
    promoteWeakToStrong(hysteresis,r+1, c+1, visited);

    return hysteresis;
}

int **createOutputPPM(int **hysteresis, string filename, int max){
       ofstream fp1(filename);
    

    if (fp1.fail()) throw("Can't open output file");

    //width and height from original input ppm file
    fp1 << "P3\n" << PPM_WIDTH << " " << PPM_HEIGHT << "\n1\n";

    for (int i = 0; i < PPM_HEIGHT ; i++)
        {
            for (int j = 0; j < PPM_WIDTH; j++) {

                if(hysteresis[i][j] == max)
                        fp1 << 1 << " " << 1 << " " << 1 << " ";
                else 
                        fp1 << 0 << " " << 0 << " " << 0 << " ";
                }
        }

        fp1.close();
        return hysteresis;
}

void createOutputPPM(const string &filename, vector<Color> &coin, int max = 255){
    ofstream out(filename);
    out << "P3" << endl << PPM_WIDTH << " " << PPM_HEIGHT << endl << max << endl;
    for(Color &d: coin){
        out << d.getR() << " " << d.getG() << " " << d.getB() << " ";
    }
    out.close();
}


void bresenLine(int **voteArr, int x1, int y1, int x2, int y2)
{
    //cout << endl << "bresenLine Begin"; 
    int x, y, dx, dy, px, py, xe, ye, i;
    dx = x2 - x1;
    dy = y2 - y1;
    int dx1 = abs(dx);
    int dy1 = abs(dy);
    px = 2 * dy1 - dx1;
    py = 2 * dx1 - dy1;
    if (dy1 <= dx1){
        if (dx >= 0){
            x = x1;
            y = y1;
            xe = x2;
        }
        else{
            x = x2;
            y = y2;
            xe = x1;
        }
        //drawPoint(x, y);
        if(x >1 && y > 1 && x < PPM_WIDTH - 1 && y < PPM_HEIGHT - 1)
        {
           //cout << endl<< "1111::" << voteArr[x][y];
        /*voteArr[x][y] = voteArr[x][y] + 1;
        if(voteArr[x][y] > MAX_VOTE)
                MAX_VOTE = voteArr[x][y];*/
        }
        
        while(x < xe)
        {
            x += 1;
            if (px < 0)
                px = px + 2 * dy1;
            else
            {
                if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
                    y += 1;
                else
                    y -= 1;
                px = px + 2 * dy1 - 2 * dx1;
            }
            //drawPoint(x, y);
           if(x >1 && y > 1 && x < PPM_WIDTH - 1 && y < PPM_HEIGHT - 1)
            {
                //cout << endl<< "2222::" << voteArr[x][y];
            /*voteArr[x][y] = voteArr[x][y] + 1;
            if(voteArr[x][y] > MAX_VOTE)
                    MAX_VOTE = voteArr[x][y];*/
            }
        }
    }
    else
    {
        if (dy >= 0)
        {
            x = x1;
            y = y1;
            ye = y2;
        }
        else
        {
            x = x2;
            y = y2;
            ye = y1;
        }
        //drawPoint(x, y);
        if(x >1 && y > 1 && x < PPM_WIDTH - 1 && y < PPM_HEIGHT - 1)
            {
               // cout << endl<< "3333::" << voteArr[x][y];
            /*voteArr[x][y] = voteArr[x][y] + 1;
            if(voteArr[x][y] > MAX_VOTE)
                    MAX_VOTE = voteArr[x][y];*/
            }
        while(y < ye)
        {
            y += 1;
            if (py <= 0)
                py = py + 2 * dx1;
            else
            {
                if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
                    x += 1;
                else
                    x -= 1;
                py = py + 2 * dx1 - 2 * dy1;
            }
            //drawPoint(x, y);
            if(x >1 && y > 1 && x < PPM_WIDTH - 1 && y < PPM_HEIGHT - 1)
            {
               // cout << endl<<  "4444::" << voteArr[x][y];
            /*voteArr[x][y] = voteArr[x][y] + 1;
            if(voteArr[x][y] > MAX_VOTE)
                    MAX_VOTE = voteArr[x][y];*/
            }
        }
    }
    //cout << endl << "bresenLine End"; 
}

Line findLine(int px, int py, double theeta, int minX, int minY, int maxX, int maxY){
    Line line;
    //radians = ( degrees * pi ) / 180 ;
    int myX1, myY1, myX2, myY2;
    bool p1Found = false;
    double thRad = (theeta * M_PI ) / 180;
    double m = tan(thRad);
    int ya = m*minX - m*px + py;
    int yb = m*maxX - m*px + py;
    int xa = (minY + m*px - py) / m;
    int xb = (maxY + m*px - py) / m;

    //cout << endl << "(0,ya):" << ya;
    //cout << endl << "(800,yb):" << yb;

    //cout << endl << "(xa,0):" << xa;
    //cout << endl << "(xb,800):" << xb;

    if(ya >=0 && ya <=800)
    {
        myY1 = ya;
        myX1 = minX;
        line.setPoint1(myX1,myY1);
         //cout << endl << "if 1::" << myX1 << "::" << myY1;
        p1Found = true;
    }
    

    if(yb >=0 && yb <=800 && p1Found == false)
    {
        myY1 = yb;
        myX1 = maxX;
        line.setPoint1(myX1,myY1);
        //cout << endl << "if 2::" << myX1 << "::" << myY1;
        p1Found = true;
    }
    else if(yb >=0 && yb <=800)
    {
        myY2 = yb;
        myX2 = maxX;
        line.setPoint2(myX2,myY2);
        //cout << endl << "if 22::" ;
        return line;     //return line
    }

    if(xa >=0 && xa <=800 && p1Found == false)
    {
        myY1 = minY;
        myX1 = xa;
        line.setPoint1(myX1,myY1);
        //cout << endl << "if 3::" << myX1 << "::" << myY1;
        p1Found = true;
    }
    else if(xa >=0 && xa <=800)
    {
        myY2 = minY;
        myX2 = xa;
        line.setPoint2(myX2,myY2);
        //cout << endl << "if 33::" ;
        return line;     //return line
    }
    
    if(xb >=0 && xb <=800)
    {
        myY1 = maxY;
        myX1 = xb;
        line.setPoint2(myX1,myY1);
        //cout << endl << "if 4::" << myX1 << "::" << myY1;
    }

    return line;

}


void incrementVote(vector<vector<int>> &vote, int x, int y, double theta)
{
    //cout << endl << "increment vote";
   /* Line l = findLine(px, py, theta, 1, 1, PPM_WIDTH - 1, PPM_HEIGHT - 1);
    cout << endl << "My point::" <<  px << "::" << py;
    if(px == 238 && py == 823)
    {
        cout << endl <<"Here it is" ;
    }
    if(&l != NULL)
        cout << endl << "Line points::" << l.getP1X()<< "::" << l.getP1Y()<<  "::" << l.getP2X() <<  "::" << l.getP2Y();
    if(l.getP1X() == 0 && l.getP1Y() == 0 && l.getP2X()  && l.getP2Y())
    {
        cout  << "All zero returning";
        return;
    }
    Point p1 = Point(l.getP1X(), l.getP1Y());
    Point p2 = Point(l.getP2X(), l.getP2Y());*/

    Point pStart;
    Point pEnd;
    double m = -tan(theta);
    if(abs(m) < 0.0001){
        pStart = Point(0, y);
        pEnd = Point(PPM_WIDTH, y);
    }else{
        double xTemp;
        pStart = Point((y/m)+x, 0);
        pEnd = Point(((-PPM_HEIGHT+y)/m)+x, PPM_HEIGHT);
    }

    int p1 = 0, q1 = 0, p2 = 0, q2 = 0, slope = 0, p = 0, q = 0, r = 0, deltaX = 0, deltaY = 0;   

    if((pEnd.getX() - pStart.getX()) > 0){
        p1 = 1;
        p2 = 1;
    }
    else{
        p1 = -1;
        p2 = -1;
    }
   if((pEnd.getY() - pStart.getY()) > 0)
        q1 = 1;
    else
        q1 = -1;

    deltaX = pEnd.getX() - pStart.getX();
    deltaY = pEnd.getY() - pStart.getY();

   p = abs(deltaX);
   q = abs(deltaY);

   //double e22 = q/p;      

     if(abs(deltaY) > abs(deltaX)){ 
        //r = abs(deltaX);
        p = abs(deltaY);
        q = abs(deltaX);
        if(deltaY > 0)
            q2 = 1;
        else
            q2 = -1;
        p2 = 0;
    }   
   slope = q/p;

   int startX = pStart.getX(); 
   int startY = pStart.getY();   
    for(int k=0;k<=p;k++){
        if(startX < PPM_WIDTH -1  && startX >=1){
            if(startY < PPM_HEIGHT -1  && startY >= 1){
            vote[startX][startY] += 1;
            if(vote[startX][startY] > MAX_VOTE)
                MAX_VOTE = vote[startX][startY];
            }
        }
        slope += q;
        if(slope>p){
            slope -= p;
            startX += p1;
            startY += q1;
        }else{
            startX += p2;
            startY += q2;
        }
    } 
    
    return;
}

void generateVoting(int **combined, vector<vector<int>> &vote, double **thetaArr)
{
    //int combinedCount = 0;
    for (int i = 1; i < PPM_HEIGHT -1 ; i++)    //**
    {
        for (int j = 1; j < PPM_WIDTH -1; j++) {    //**
            if(combined[i][j] == 1)     //for each edge
            {
                //combinedCount++;
                incrementVote(vote, i, j, thetaArr[i][j]);
            }
        }
    }
    //cout <<endl << "combinedCount:" << combinedCount;
}

void createVotePPM(string &filename, vector<vector<int>> &vote, int vMax = 1){
    ofstream out(filename);
    out << "P3" << endl << PPM_WIDTH << " " << PPM_HEIGHT << endl << vMax << endl;
    int m;
    for(size_t i=0;i<vote[0].size();i++){
        for(size_t j=0;j<vote.size();j++){
            if( vote[j][i] > vMax)
                out << vMax << " " << vMax << " " << vMax << " ";
            else
                out << vote[j][i] << " " << vote[j][i] << " " << vote[j][i] << " ";
        }
    }
    out.close();
}

void createVotePPM1(string &filename, int **vote, int vMax = 1){
    ofstream out(filename);
    out << "P3" << endl << PPM_WIDTH << " " << PPM_HEIGHT << endl << vMax << endl;
    int m;
    for(size_t i=0;i<PPM_HEIGHT ;i++){
        for(size_t j=0;j<PPM_WIDTH;j++){
            if( vote[i][j] > vMax)
                out << vMax << " " << vMax << " " << vMax << " ";
            else
                out << vote[i][j] << " " << vote[i][j] << " " << vote[i][j] << " ";
        }
    }
    out.close();
}

void createPPMs(int **arr, vector<Color> *imputPPM)
{
    int** a = new int*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        a[i] = new int[PPM_WIDTH];

    int** hysteresis = new int*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        hysteresis[i] = new int[PPM_WIDTH];

    double** g_x = new double*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        g_x[i] = new double[PPM_WIDTH];
    
    double** g_y = new double*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        g_y[i] = new double[PPM_WIDTH];

    int** nms = new int*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        nms[i] = new int[PPM_WIDTH];

    int** combined = new int*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        combined[i] = new int[PPM_WIDTH];

    double** thetaArr = new double*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        thetaArr[i] = new double[PPM_WIDTH];

    int** voteArr = new int*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        voteArr[i] = new int[PPM_WIDTH];

    vector<vector<int>> vote(PPM_WIDTH, vector<int>(PPM_HEIGHT, 0));

    string filenameHyst = "image2.ppm"; 
    string filenameNms = "image1.ppm"; 
    string filenameCombine = "imagef.ppm"; 
    string filenameVote = "imagev.ppm"; 
    string filenameCoin = "coins.ppm";

    double gX = 0.0, gY = 0.0, g = 0.0;
    double theta = 0.0;

    for (int i = 1; i < PPM_HEIGHT -1 ; i++)
    {
        for (int j = 1; j < PPM_WIDTH -1; j++) {
            //try1 - column major
            gX = 0.0;
            gY = 0.0;
            //gX = ( (-1 * arr[j-1][i-1]) + (arr[j+1][i-1]) + (-2*arr[j-1][i]) + (2*arr[j+1][i]) + (-1*arr[j-1][i+1]) + (arr[j+1][i+1]) );
            //gY = ( (-1 * arr[j-1][i-1]) + (-2 * arr[j][i-1]) + (-1*arr[j+1][i-1]) + (arr[j-1][i+1]) + (2*arr[j][i+1]) + (arr[j+1][i+1]) );
            //try2 - row major - works
            gX = ( (-1 * arr[i-1][j-1]) + (arr[i+1][j-1]) + (-2*arr[i-1][j]) + (2*arr[i+1][j]) + (-1*arr[i-1][j+1]) + (arr[i+1][j+1]) );
            gY = ( (-1 * arr[i-1][j-1]) + (-2 * arr[i][j-1]) + (-1*arr[i+1][j-1]) + (arr[i-1][j+1]) + (2*arr[i][j+1]) + (arr[i+1][j+1]) );

            g_x[i][j] = gX;
            g_y[i][j] = gY;

            g = sqrt( (gX * gX) + (gY * gY));
            theta = atan2(gY, gX);
            thetaArr[i][j] = theta;

           /*if(g < LOW_THRES)
                hysteresis[i][j] = 0;
            else if(g >= LOW_THRES && g <= HIGH_THRES)
                hysteresis[i][j] = 1;
            else
                hysteresis[i][j] = 2;*/

            if(g > HIGH_THRES)
                hysteresis[i][j] = 2;
            else if(g >= LOW_THRES)
                hysteresis[i][j] = 1;
            else
                hysteresis[i][j] = 0;

            nms[i][j] = 0; //initialize with 0
            voteArr[i][j] = 0; //initialize with 0

            if ( ( (theta < 22.5) && (theta > -22.5) ) || (theta > 157.5) || (theta < -157.5) ) //0 degree
            {
				if(arr[i][j] >= arr[i][j-1] && arr[i][j] >= arr[i][j+1])
                    nms[i][j] = 1;
            }
			if ( ( (theta > 22.5) && (theta < 67.5) ) || ( (theta < -112.5) && (theta > -157.5) ) ) //45 degree
            {
				if(arr[i][j] >= arr[i-1][j+1] && arr[i][j] >= arr[i+1][j-1])
                    nms[i][j] = 1;
            }
			if ( ( (theta > 67.5) && (theta < 112.5) ) || ( (theta < -67.5) && (theta > -112.5) ) ) //90 degree
            {
				if(arr[i][j] >= arr[i+1][j] && arr[i][j] >= arr[i-1][j])
                    nms[i][j] = 1;
            }
			if ( ( (theta > 112.5) && (theta < 157.5) ) || ( (theta < -22.5) && (theta > -67.5) ) ) //135 degree
            {
				if(arr[i][j] >= arr[i+1][j+1] && arr[i][j] >= arr[i-1][j-1])
                    nms[i][j] = 1;
            }

        }

    }

    bool** visited = new bool*[PPM_HEIGHT];
    for(int i = 0; i < PPM_HEIGHT; ++i)
        visited[i] = new bool[PPM_WIDTH];

    for (int i = 0; i < PPM_HEIGHT ; i++)
    {
        for (int j = 0; j < PPM_WIDTH; j++) {
              visited[i][j] = false;  
        }
    }

    for (int i = 1; i < PPM_HEIGHT -1 ; i++)
    {
        for (int j = 1; j < PPM_WIDTH -1; j++) {
            if(hysteresis[i][j] == 2)
               hysteresis = promoteWeakToStrong(hysteresis, i, j, visited);
        }
    }

    for (int i = 1; i < PPM_HEIGHT -1 ; i++)
    {
        for (int j = 1; j < PPM_WIDTH -1; j++) {
            if(hysteresis[i][j] == 2 && nms[i][j] == 1)
               combined[i][j] = 1;
            else 
                combined[i][j] = 0;
        }
    }

    createOutputPPM(hysteresis, filenameHyst, 2);
    createOutputPPM(nms, filenameNms, 1);

    //merge Hysteresis + NMS
    for (int i = 1; i < PPM_HEIGHT -1 ; i++)
    {
        for (int j = 1; j < PPM_WIDTH -1; j++) {
            combined[i][j] = 0; //initialize
            if(hysteresis[i][j] == 2 && nms[i][j] == 1)
               combined[i][j] = 1;
        }
    }

    createOutputPPM(combined, filenameCombine, 1);

   // cout << "BEFORE MAX_VOTE: " << MAX_VOTE;

    ////generateVoting(combined, vote, thetaArr);

    //cout << "AFTER MAX_VOTE: " << MAX_VOTE;

    //createOutputPPM(voteArr, filenameVote, MAX_VOTE);

    ////createVotePPM(filenameVote, vote, MAX_VOTE);

    //cout << "After File MAX_VOTE: " << MAX_VOTE;

    /*voteArr[i][j] = voteArr[i][i] + 1;
        if(voteArr[i][j] > MAX_VOTE)
                MAX_VOTE = voteArr[i][j];*/

       // Circle Detection
    int rmin = 5;
    int rmax = 60;
    // double theta_step = 5.0;
    int circ_threshold = 6;

    int curr_index;
    int a_, b_;
    double t;
    map<string, int> acc = {};
    string st;
    for (int x = 1; x < PPM_HEIGHT-1; x++) {
        for (int y = 1; y < PPM_WIDTH-1; y++) {
            curr_index = PPM_WIDTH*x+y;
            if (combined[x][y] == 1) {
                for (int sign = -2; sign < 2; sign+=1) {
                    t = atan2(g_y[x][y], g_x[x][y]) + (double) sign * M_PI / 2.0;
                    for (int r = rmin; r <= rmax; r++) {
                        a_ = x - (int) r * cos(t);
                        b_ = y - (int) r * sin(t);
                        if (a_ < 0 || a_ >= PPM_HEIGHT-1 || b_ < 0 || b_ >= PPM_WIDTH-1) break;
                        st = to_string(a_) + " " + to_string(b_) + " " + to_string(r);
                        voteArr[a_][b_] = voteArr[a_][b_] + 1;
                        if(voteArr[a_][b_] > MAX_VOTE)
                            MAX_VOTE = voteArr[a_][b_];
                        acc[st]++;
                    }
                }
            }
        }
    }
    //cout << endl << "voting size: " << acc.size();

    //cout << endl << "MAX VOTE: " << MAX_VOTE;

    createVotePPM1(filenameVote, voteArr, MAX_VOTE);

    vector<pair<string, int>> vecVotes;
    for (auto i : acc) vecVotes.push_back(make_pair(i.first, i.second));
    sort(vecVotes.begin(), vecVotes.end(), sortbysecdesc);

    unordered_set<int> cirPos;
    map<int, int> radPos;
    map<int, int> coinType;
    int x, y, ch_r;
    int itr;
    int r = 0;
    vector<int> abr;
    int frq = 0;
    bool circle_check;
    int min_rad = max(PPM_WIDTH, PPM_HEIGHT);
    for (auto i : vecVotes) {
        st = i.first;
        frq = i.second;
        abr.clear();
        istringstream ss(st);
        for (itr = 0; itr < 3; itr++) {
            string word;
            ss >> word;
            abr.push_back(stoi(word));
            // cout << word << " ";
        }
        // cout << endl;
        a_ = abr[0]; b_ = abr[1]; r = abr[2];
        if (frq >= circ_threshold) {
            circle_check = true;
            for (auto ch : radPos) {
                x = ch.first / PPM_WIDTH;
                y = ch.first % PPM_WIDTH;
                ch_r = ch.second;
                if ((a_ - x)*(a_ - x) + (b_ - y)*(b_ - y) <= ch_r * ch_r) {
                    circle_check = false;
                        break;
                }
            }
            if (!circle_check) continue;
            // cout << frq << " " << a_ << " " << b_ << " " << r << endl;
            cirPos.insert(PPM_WIDTH * a_ + b_);
            radPos[PPM_WIDTH * a_ + b_] = r;
            min_rad = min(min_rad, r);
        }
    }

    //cout << endl << " Total centers: " << cirPos.size();
    //cout << endl << " Total radPos: " << radPos.size();

   /* int index;
    for(int i=0;i<PPM_HEIGHT-1;i++){
        for(int j=0;j<PPM_WIDTH-1;j++){
            index = PPM_WIDTH*i+j;
            if (cirPos.find(index) != cirPos.end()){
               // cout << endl << "found center";
                imputPPM->at(index) = Color(255, 0, 0);
                for(int r=0;r<5;r++)
                    createCircle(*imputPPM, i, j, r, 255, 0, 0);
            }
        }
    } */

    int cents = 0;
    vector<double> diam_ratios =
            {750. / 705., 835. / 705., 705. / 705., 955. / 705., 1205. / 705., 1043. / 705.}; // p, n, d, q, hD, D

    min_rad = 27;
    int penn = 0; int nick = 0; int dime = 0; int quar = 0; int half = 0; int doll = 0;
    int index = 0; double min_diff; double this_diff;
    for (auto i : radPos) {
        ch_r = i.second;
        min_diff = (double) max(PPM_WIDTH, PPM_HEIGHT);
        for (itr = 0; itr < 6; itr++) {
            this_diff = abs(ch_r / (double) min_rad - diam_ratios[itr]);
            //cout << endl << "ch_r: " << ch_r << " min_rad: " << min_rad << " diam_ratios[itr]: " << diam_ratios[itr] << " this_diff: " << this_diff;
            if (this_diff < min_diff) {
                min_diff = this_diff;
                index = itr;
                //cout << endl << "index:" << index;
            }
        }
        //cout << endl << " index final: " << index;
        if (index == 0) 
            { penn++;  
                coinType.insert(pair<int, int>(i.first, 1));
            }
        else if (index == 1) 
        {
            nick++;  
            coinType.insert(pair<int, int>(i.first, 2));

        }
        else if (index == 2) 
        {
            dime++;  
            coinType.insert(pair<int, int>(i.first, 3));
        }
        else if (index == 3) 
        {
                quar++;
                coinType.insert(pair<int, int>(i.first, 4));
        }
        else if (index == 4) 
        {
                half++;
                coinType.insert(pair<int, int>(i.first, 5));
        }
        else if (index == 5)
        {
                doll++;
                coinType.insert(pair<int, int>(i.first, 6));
        }
    }

    coinType.size();
     int k;
    for(int i=0;i<PPM_WIDTH;i++){
        for(int j=0;j<PPM_HEIGHT;j++){
            k = j*PPM_WIDTH+i;
            if (cirPos.find(k) != cirPos.end()){
            imputPPM->at(k) = Color(255, 0, 0);
                for(int r=0;r<5;r++)
                    createCircle(*imputPPM, i, j, r, 255, 0, 0);
                 map<int, int>::iterator itr;
                 itr = coinType.find(k);
                    if(itr->second == 1)
                    {
                        //cout << "Penny stored:" << endl;
                        for(int r=0;r<5;r++){
                            int coinR = (int) (28 + r);
                            createCircle(*imputPPM, i, j, coinR, 255, 0, 0);
                        }
                    }
                    if(itr->second == 2)
                    {
                        //cout << "nick stored:" << endl;
                        for(int r=0;r<5;r++){
                            int coinR = (int) (38 + r);
                            createCircle(*imputPPM, i, j, coinR, 0, 0, 255);
                        }
                    }
                    if(itr->second == 3)
                    {
                        //cout << "dime stored:" << endl;
                        for(int r=0;r<5;r++){
                            int coinR = (int) (24 + r);
                            createCircle(*imputPPM, i, j, coinR, 0, 255, 0);
                        }
                    }
                    if(itr->second == 4)
                    {
                        //cout << "qtr stored:" << endl;
                        for(int r=0;r<5;r++){
                            int coinR = (int) (40 + r);
                            createCircle(*imputPPM, i, j, coinR, 255, 255, 0);
                        }
                    }
                    if(itr->second == 5)
                    {
                        //cout << "half stored:" << endl;
                        for(int r=0;r<5;r++){
                            int coinR = (int) (45 + r);
                            createCircle(*imputPPM, i, j, coinR, 0, 255, 255);
                        }
                    }
                    if(itr->second == 6)
                    {
                        //cout << "dollar stored:" << endl;
                        for(int r=0;r<5;r++){
                            int coinR = (int) (45 + r);
                            createCircle(*imputPPM, i, j, coinR, 0, 255, 255);
                        }
                    }

            }
        }
    }

    int sum = penn + 5 * nick + 10 * dime + 25 * quar + 50 * half + 100 * doll;
    cout << endl << "Total Coins:: " << "penny:" << penn << " dime:" << dime << " quar:" << quar << " half:" << half << " doll:" << doll;
    cout << endl << "Total: " << sum << " cents" << endl;

    ofstream out_file;
	out_file.open("results.txt");
	out_file << "Total Coins:: " << "penny:" << penn << " dime:" << dime << " quar:" << quar << " half:" << half << " doll:" << doll;
    out_file << endl << "Total: " << sum << " cents" << endl;

	out_file.close(); 

    createOutputPPM(filenameCoin, *imputPPM);

    return;
}

//From vector result to mapped kmeans to create output ppm
void createPPMOutput(int **arr)
{
    string filename = "imageg.ppm"; 
    ofstream fp(filename);

    if (fp.fail()) throw("Can't open output file");

    //width and height from original input ppm file
    fp << "P3\n" << PPM_WIDTH << " " << PPM_HEIGHT << "\n255\n";

    for (int i = 0; i < PPM_HEIGHT; i++)
    {
        for (int j = 0; j < PPM_WIDTH; j++) {

             fp <<  arr[i][j] << " " <<arr[i][j] << " " << arr[i][j] << " ";

        }
        /*if(count > 81)
                break;*/
    }
    fp.close();
}



void readPpmFile()
{
    int width, height;
    string filename = "image.ppm"; 
    vector<Color> *inputPPM = new vector<Color>();
    ifstream fp(filename);
    if (fp.fail())
    {
        cout << "inout file not found";
        return; //You failed
    }

    //Read the Magic Number
    string mg_num, width_str, height_str, range_str;
    fp >> mg_num;

    if (mg_num != "P3") {
        fp.close();
        return; //The file is not a ASCII PPM file!
    }

    fp >> width_str >> height_str >> range_str;
    width  = atoi(width_str .c_str()); //Takes the number string and converts it to an integer
    height = atoi(height_str.c_str());
    PPM_WIDTH = width;
    PPM_HEIGHT = height;

    //cout << "WIDTH:HEIGHT:" << PPM_WIDTH << ":" << PPM_HEIGHT;

    //RGB tmp;
    Color tmp;
    string _R, _G, _B;

    int** ppmArr = new int*[height];
    for(int i = 0; i < height; ++i)
        ppmArr[i] = new int[width];

    //int count = 0;

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++) {
            //fin >> colors.r >> colors.g >> colors.b;
            fp >> _R >> _G >> _B;
            tmp.setColor(atoi(_R.c_str()), atoi(_G.c_str()), atoi(_B.c_str()));
            inputPPM->push_back(Color(tmp.getR(), tmp.getG(), tmp.getB()));
            int avg = (tmp.getR() + tmp.getG() + tmp.getB() ) /3;
            //cout << "(" << i << ", " << j << "): " << colors.r << " " << colors.g << " " << colors.b << endl;

            //Column major - row index comes first.
            ppmArr[i][j] = avg;   //Copy RGB values into image
           // count++;
        }
    }
   // cout << "count" << count;

    fp.close();

    //createPPMOutput(ppmArr);    //generate gray image ppm

    //createBinaryPPM(ppmArr);    //generate binary image ppm
    createPPMs(ppmArr, inputPPM);

    
}

//From vector result to mapped kmeans to create output ppm
void createPPMOutput(std::vector<int> *vect_r, std::vector<Color> *vect_m)
{
    string filename = "output.ppm"; 
    ofstream fp(filename);

    if (fp.fail()) throw("Can't open output file");

    //width and height from original input ppm file
    fp << "P3\n" << PPM_WIDTH << " " << PPM_HEIGHT << "\n255\n";

    for (auto it = vect_r->begin(); it != vect_r->end(); ++it) 
    {
        int ind = *it;
        //cout << endl << ind << "(" <<  vect_m->at(ind).getR() << "," << vect_m->at(ind).getG() << "," << vect_m->at(ind).getB() << ")";
        fp <<  vect_m->at(ind).getR() << " " << vect_m->at(ind).getG() << " " << vect_m->at(ind).getB() << " ";
        //cout << endl << std::fixed << std::setprecision(23) << "(" << c.getR() << "," << c.getG() << "," << c.getB() << ")"  ;
    }
    //cout << "Result count: " << vect_r->size();
}

void readFromFile(std::vector<Point> *vect)
{
    // filestream variable file 
    fstream file; 
    string filename; 
    char cPoint[200];
    string sPoint;
    string* points = new string[10];
    int i = 0;
  
    // filename of the file 
    filename = "points100.txt"; 
  
    // opening file 
    file.open(filename.c_str()); 
    if(!file) {
        cout << "Cannot open input file points.txt.\n";
        return;
    }
  
    // extracting words from the file 
    while (file.good()) 
    { 
        // displaying content 
        string intermediate; 
        file.getline(cPoint, 256); 
        //int pSize = sizeof(cPoint) / sizeof(char); 
        sPoint = convertToString(cPoint); 
        stringstream sLine(sPoint); 
        
        double x, y;
        int odd = 1;
        sPoint.erase(remove(sPoint.begin(), sPoint.end(), ' '), sPoint.end());
        if(!sPoint.empty())
        {
            while(getline(sLine, intermediate, ' '))
            { 
                intermediate.erase(remove(intermediate.begin(), intermediate.end(), ' '), intermediate.end());
                if(!intermediate.empty())
                {
                    if(odd % 2 != 0)
                        x = stod(intermediate);
                    else
                        y = stod(intermediate);
                    
                    points[i] = intermediate;
                    odd++;
                }
            } 
            Point p(x,y);
            //cout << endl << std::setprecision(23) << "point" << p.getX() << ":" <<  p.getY(); 
            vect->push_back(p);
        }
    } 
    file.close();
    return;
}

std::vector<Point> createPontsList(list<double> points)
{
    std::vector<Point> vec_points;
    for (auto it = points.begin(); it != points.end(); ++it) 
    {
        double x = *it;
        ++it;
        double y = *it;
        Point pt = Point(x, y);
        vec_points.push_back(pt);
        //cout << endl << std::fixed << std::setprecision(23) << pt.getX() << "::" << pt.getY();
    }
    return vec_points;
}

std::vector<Color> generateColorList(int count)
{
    std::vector<Color> vec_points;

    for(int i = 0; i < count ; i++)
    {
        vec_points.push_back(col[i]);
    }
    return vec_points;
}

std::vector<Color> initializeMeans(int count)
{
    std::vector<Color> vec_points;

    for(int i = 0; i < count ; i++)
    {
        vec_points.push_back(col[i]);
    }
    return vec_points;
}

/*void createCircles(std::vector<Point> *vect, std::vector<Point> *vect1, std::vector<Color> *vect_col, std::vector<int> *vect_r)
{
    int result_index = 0;
    for (auto it1 = vect->begin(); it1 != vect->end(); ++it1) 
    {
        Point p = *it1;
        int cluster = vect_r->at(result_index);
        Color c = vect_col->at(cluster);
        createCircle(p.getX()*800.0, p.getY()*800.0, 2,c.getR(), c.getG(), c.getB());
        result_index++;
    }

    int i = 0;
    for (auto it1 = vect1->begin(); it1 != vect1->end(); ++it1) 
    {
        Point p = *it1;
        //Color col = vect_col->at(i);
        createCircle(p.getX()*800.0, p.getY()*800.0, 3, 0, 0, 0);
        i++;
    }
}*/

void part2()
{
    srand(time(nullptr));
    readPpmFile();
}

int main() { 
    part2();
}