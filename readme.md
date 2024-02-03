# Fixed-wing Simulator

用于固定翼无人机的轨迹规划模块仿真

* 路径搜索 A*
* 轨迹生成 Bspling
* 轨迹重优化
* 由传感器范围有限所导致的重规划

## 安装

安装系统依赖
```
sudo apt-get install cmake libopenblas-dev liblapack-dev libarpack-dev libarpack2-dev libsuperlu-dev
```

安装Armadillo
```
xz -d armadillo-9.870.2.tar.xz
tar -xvf armadillo-9.870.2.tar
cd armadillo-9.870.2
mkdir build
cd build
cmake ..
make
sudo make install
```
## 功能包介绍

* random_complex：随机生成障碍物点云地图；
* waypoint_generator：给定目标点；
* odom_visualization：四旋翼可视化；
* pcl_render_node：简单版的局部传感器模型，返回局部范围内的障碍物点云；
* **trajectory_generator_node** ：规划一条可行的Bspline轨迹；
* traj_server：将Bspline轨迹转换为控制指令；
* so3_control：将控制指令转换为实际控制量；
* quadrotor_simulator_so3：固定翼无人机仿真模型。

