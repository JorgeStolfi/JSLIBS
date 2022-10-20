/* See {}. */

/* Exports planar slices */

#include <stdint.h>

int32_t exportSingleSVGFormat (string fileName, const vector<vector<salamic_r2_Segment_t> > &slicesWithLineSegs, const salamic_stl_r3_t &aabbSize) {

    typedef struct _rgb {
       int32_t r, g, b;
    } rgb; 
  
    int32_t nrgb = 8;   
 
    rgb colors[nrgb];

    colors[0] = {128,   0, 128}; /*Purple*/
    colors[1] = {  0,   0,   0}; /*Black*/
    colors[2] = {255,   0,   0}; /*Red*/
    colors[3] = {  0, 128,   0}; /*Green*/
    colors[4] = {  0,   0, 255}; /*Blue*/
    colors[5] = {  0, 255, 255}; /*Cyan*/
    colors[6] = {128, 128,   0}; /*Olive*/
    colors[7] = {128, 128, 128}; /*Gray*/

    //char filename[256];
    FILE *f=NULL;
    float dx=0, dy=0;
    f=fopen(fileName.c_str(), "w");
	//fprintf(f,"<html>\n");
	//fprintf(f,"<body>\n");
	fprintf(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	fprintf(f, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
	fprintf(f, "<svg preserveAspectRatio=\"xMidYMid meet\" width=\"1024\" height=\"768\" viewBox=\"0 0 1024 768\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:cc=\"http://web.resource.org/cc/\">\n");
    if (!f)    return 1;
    const size_t nSlices = slicesWithLineSegs.size();
    const size_t slicePerRow = (size_t)sqrt((float)nSlices)*2;
    for (size_t i=0; i<nSlices; ++i) {
		
     	dx = (float)(i%slicePerRow)*(aabbSize.x*1.05f);
        dy = (float)(i/slicePerRow)*(aabbSize.y*1.05f);

        for (const salamic_r2_Segment_t &ls : slicesWithLineSegs[i]) {

           fprintf(f, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(%d,%d,%d)\"/>\n",
                          dx + ls.v[0].x, 
                          dy + ls.v[0].y,
                          dx + ls.v[1].x, 
                          dy + ls.v[1].y, 
                          colors[ls.label % nrgb].r, 
                          colors[ls.label % nrgb].g, 
                          colors[ls.label % nrgb].b);
 	}
    }
    fprintf(f,"</svg>\n");
    fclose(f);
    return 0;
}

void exportSingleSVGFormatB (const vector<vector<salamic_r2_Segment_t> > &slicesWithLineSegs, const salamic_stl_r3_t &aabbSize, int32_t nsegments, int32_t nslices) {

    glm::vec3 fromEuler (0.0f, 0.0f, 60.0f);
    glm::quat quaternion (DEG_TO_RAD(fromEuler));
    glm::vec3 toEuler = glm::eulerAngles(quaternion);
    float angle = glm::angle(quaternion);
    glm::vec3 axis = glm::axis(quaternion);
    glm::mat4 View = glm::rotate(glm::mat4(1.0), angle, axis);
    float zoom = 0.3f;
    glm::mat4 Projection = glm::perspective (zoom, 4.0f / 3.0f, 0.1f, 100.f);
    glm::mat4 Model = glm::lookAt (
        glm::vec3(1, 1, 1),    // Eye point (where am I?)
        glm::vec3(0, 0, 0),    // Look at Center point (What am I looking at?)
        glm::vec3(0, 1, 0));   // UP vector of the camera (Where is my UP vector?)
    glm::mat4 MVP = Projection * View * Model;

    // Start Output
    char filename[256];

    sprintf(filename, "./video/%05d_%05d.svg", nslices, nsegments);

    FILE *file = fopen(filename, "w");

    fprintf (file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf (file, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
    fprintf (file, "<svg viewBox=\"0 0 1024 768\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:cc=\"http://web.resource.org/cc/\">\n");

    size_t i = 0;
    int32_t disp_x = +100; /*Deslocamento da figura*/
    int32_t disp_y = +100; /*Deslocamento da figura*/
    //int32_t shift_x = -180; /*Bottle*/
    //int32_t shift_y = -100; /*Bottle*/
    int32_t shift_x = +150;
    int32_t shift_y = +250; /*Quanto mais positivo, mais abaixa a camera!*/

    for (; i < nslices; i++) {
       for (const salamic_r2_Segment_t &ls : slicesWithLineSegs[i]) {
          salamic_stl_r3_t p0 = disp_x + ls.v[0];
          salamic_stl_r3_t p1 = disp_y + ls.v[1];
          p0.transform(MVP);
          p1.transform(MVP);
          fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(0,0,1)\"/>\n", p0.x + shift_x, p0.y + shift_y, p1.x + shift_x,p1.y + shift_y);
       }
    }

    int32_t segment = 0;
    for (const salamic_r2_Segment_t &ls : slicesWithLineSegs[i]) {
       if (segment < nsegments) {
           salamic_stl_r3_t p0 = disp_x + ls.v[0];
           salamic_stl_r3_t p1 = disp_y + ls.v[1];
           p0.transform(MVP);
           p1.transform(MVP);
           fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(0,0,1)\"/>\n", p0.x + shift_x, p0.y + shift_y, p1.x + shift_x, p1.y + shift_y);
           segment++;
        }
    } 

    fprintf (file,"</svg>\n");

    fclose(file);
}

int32_t exportSingleSVGFormat3D (const vector<vector<salamic_r2_Segment_t> > &slicesWithLineSegs, const salamic_stl_r3_t &aabbSize) {

    const size_t nslices = slicesWithLineSegs.size();

    //printf("Número de fatias: %d\n", nslices);
    for (size_t i = 0; i < nslices; i++) {
        const size_t nsegments = slicesWithLineSegs[i].size();
        for (size_t j = 0; j < nsegments; j++) {
           exportSingleSVGFormatB (slicesWithLineSegs, aabbSize, j, i);
        }
    }
    return 0;
}

int32_t exportSingleSVGFormat4D (const vector<vector<salamic_r2_Segment_t> > &slicesWithLineSegs,  const salamic_stl_r3_t &aabbSize) {

    glm::vec3 fromEuler (0.0f, 0.0f, 60.0f);
    glm::quat quaternion (DEG_TO_RAD(fromEuler));
    glm::vec3 toEuler = glm::eulerAngles(quaternion);
    float angle = glm::angle(quaternion);
    glm::vec3 axis = glm::axis(quaternion);
    glm::mat4 View = glm::rotate(glm::mat4(1.0), angle, axis);
    float zoom = 0.3f;
    glm::mat4 Projection = glm::perspective (zoom, 4.0f / 3.0f, 0.1f, 100.f);
    glm::mat4 Model = glm::lookAt (
        glm::vec3(1, 1, 1),    // Eye point (where am I?)
        glm::vec3(0, 0, 0),    // Look at Center point (What am I looking at?)
        glm::vec3(0, 1, 0));   // UP vector of the camera (Where is my UP vector?)
    glm::mat4 MVP = Projection * View * Model;


    // Start Output
    char filename[256];
    FILE *f=NULL;
    salamic_stl_r3_t p0, p1;
    f=fopen("out3D.svg", "w");
    if (!f)    return 1;
	fprintf(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	fprintf(f, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
	fprintf(f, "<svg viewBox=\"0 0 1024 768\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:cc=\"http://web.resource.org/cc/\">\n");
    const size_t nSlices = slicesWithLineSegs.size();
    const size_t slicePerRow = (size_t)sqrt((float)nSlices);
    for (size_t i=0; i<nSlices; ++i) {
        //const vector<vector<LineSegment> > &contour = slicesWithLineSegs[i];

        for (const salamic_r2_Segment_t &ls : slicesWithLineSegs[i]) {
            //fprintf(f, "    <polygon points=\"");
			p0=+100+ls.v[0];
			p1=+100+ls.v[1];
			p0.transform(MVP);
			p1.transform(MVP);
			fprintf(f, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(0,0,1)\"/>\n",
				p0.x,p0.y,p1.x,p1.y);
			//fprintf(f, "%f,%f %f,%f ", p0.x,p0.y,p1.x,p1.y);
            //fprintf(f, "\" style=\"stroke:blue;stroke-width:1\" \\>\n");
        }
    }
	fprintf(f,"</svg>\n");
    fclose(f); 

 
    //glm::vec3 angle = glm::eulerAngles (100.0f, 30.0f, 60.0f); /*Bottle*/

    //glm::detail::tvec3<float> eulerAngles(100.0f, 30.0f, 60.0f); /*Bottle*/
    //glm::detail::tvec3<float> eulerAngles(0.0f, 0.0f, 60.0f); /*Bottle*/
    //glm::detail::tvec3<float> eulerAngles(0.0f, 0.0f, 60.0f); /*Soldier*/
    //glm::detail::tvec3<float> eulerAngles(90.0f, 0.0f, 0.0f); /*Femur*/

    //glm::detail::tquat<float> quaternion = glm::detail::tquat<float>(DEG_TO_RAD(toEuler)); // This glm func wants values as radians
    /*float angle = glm::angle(quaternion);
    glm::detail::tvec3<float> axis = glm::axis(quaternion);

    glm::mat4 View        = glm::rotate(glm::mat4(1.0), angle, axis);
    glm::mat4 Projection  = glm::perspective(2.0f, 4.0f / 3.0f, 0.1f, 100.f);
    //glm::mat4 Projection  = glm::perspective(5.0f, 4.0f / 3.0f, 0.1f, 100.f);
    //glm::mat4 Projection  = glm::perspective(20.0f, 4.0f / 3.0f, 0.1f, 100.f);
    glm::mat4 Model       = glm::lookAt(    
        glm::vec3(1, 1, 1),    // Eye point (where am I?)
        glm::vec3(0, 0, 0),    // Look at Center point (What am I looking at?)
        glm::vec3(0, 1, 0));   // UP vector of the camera (Where is my UP vector?)
    glm::mat4 MVP = Projection * View * Model; 

    // Start Output
    char filename[256];
    FILE *f=NULL;
    salamic_stl_r3_t p0, p1;
    f=fopen("out3D.svg", "w");
    if (!f)    return 1;
	fprintf(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	fprintf(f, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
	fprintf(f, "<svg viewBox=\"0 0 1024 768\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:cc=\"http://web.resource.org/cc/\">\n");
    const size_t nSlices = slicesWithLineSegs.size();
    const size_t slicePerRow = (size_t)sqrt((float)nSlices);
    for (size_t i=0; i<nSlices; ++i) {
        //const vector<vector<salamic_r2_Segment_t> > &contour = slicesWithLineSegs[i];

        for (const salamic_r2_Segment_t &ls : slicesWithLineSegs[i]) {
            //fprintf(f, "    <polygon points=\"");
			p0=+100+ls.v[0];
			p1=+100+ls.v[1];
			p0.transform(MVP);
			p1.transform(MVP);
			fprintf(f, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(0,0,1)\"/>\n",
				p0.x,p0.y,p1.x,p1.y);
			//fprintf(f, "%f,%f %f,%f ", p0.x,p0.y,p1.x,p1.y);
            //fprintf(f, "\" style=\"stroke:blue;stroke-width:1\" \\>\n");
        }
    }
	fprintf(f,"</svg>\n");
    fclose(f); */
    return 0;
}

