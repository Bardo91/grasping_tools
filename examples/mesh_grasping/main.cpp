//---------------------------------------------------------------------------------------------------------------------
//  GRASPING_TOOLS
//---------------------------------------------------------------------------------------------------------------------
//  Copyright 2018 Pablo Ramon Soria (a.k.a. Bardo91) pabramsor@gmail.com
//---------------------------------------------------------------------------------------------------------------------
//  Permission is hereby granted, free of charge, to any person obtaining a copy of this software
//  and associated documentation files (the "Software"), to deal in the Software without restriction,
//  including without limitation the rights to use, copy, modify, merge, publish, distribute,
//  sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in all copies or substantial
//  portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
//  BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
//  OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//---------------------------------------------------------------------------------------------------------------------

#include <iostream>

#include <iostream>
#include <thread>
#include <chrono>

#include <pcl/visualization/pcl_visualizer.h>

#include <grasping_tools/grasping_tools.h>
#include <gpis/sampler/prior/PriorFactory.h>
#include <gpis/utils/SurfaceGpis.h>

bool newGrasp = false;
void keycallback(const pcl::visualization::KeyboardEvent &_event, void *_data) {
    if (_event.keyDown() && _event.getKeySym() == "n") {
        newGrasp = true;
    }
};

int main(int _argc, char** _argv){

    if(_argc != 3){
        std::cout << "Bad input arguments. Usage: ./example_mesh_grasping [resolution] [path to stl mesh file]" << std::endl;
        return -1;
    }

    std::cout << "Creating object" << std::endl;

    float resolution = atof(_argv[1]);
    grasping_tools::ObjectMesh object(_argv[2]);
    pcl::PolygonMesh mesh;
    object.mesh(mesh);

    std::cout << "Creating candidate points" << std::endl;
    auto centroFaces = object.centroidFaces();

    std::cout << "Creating Gripper" << std::endl;
    grasping_tools::GripperHandTwoFingers gripper(1.0);

    std::cout << "Entering in the loop" << std::endl;
    auto grasps = gripper.generateGrasps(object, resolution);
    std::cout << "Found: " << grasps.size() << " possible grasps" << std::endl;
    auto grasp = grasps.begin();
    pcl::PointCloud<pcl::PointXYZRGB> cloudUsedPoints;
    int counter = 0; 
    
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->registerKeyboardCallback(keycallback, nullptr);
    viewer->setBackgroundColor(1,1,1);
    viewer->addPolygonMesh(mesh, "object");

    for(unsigned i = 0; i < grasps.size() ; i++){
        auto cps = grasps[i].contactPoints();

        viewer->removeAllShapes();

        int counter = 0;
        for(auto &cp: cps){
            plotContactPoint(cp, *viewer, 1, "cone"+std::to_string(counter++));
        }
    
        std::cout << "Plotting grasp " << i << std::endl;
        while(!newGrasp)
            viewer->spinOnce(30,true);

        newGrasp = false;
    }
    
    viewer->spin();

}