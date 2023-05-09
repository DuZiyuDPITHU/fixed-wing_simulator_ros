#include "Astar_searcher.h"

using namespace std;
using namespace Eigen;

inline double norm2(Vector3i a, Vector3i b)
{
  return sqrt(double((a(0)-b(0))*(a(0)-b(0))+(a(1)-b(1))*(a(1)-b(1))+(a(2)-b(2))*(a(2)-b(2))));
}

void AstarPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l,
                                  Vector3d global_xyz_u, int max_x_id,
                                  int max_y_id, int max_z_id) {
  gl_xl = global_xyz_l(0);
  gl_yl = global_xyz_l(1);
  gl_zl = global_xyz_l(2);

  gl_xu = global_xyz_u(0);
  gl_yu = global_xyz_u(1);
  gl_zu = global_xyz_u(2);

  GLX_SIZE = max_x_id;
  GLY_SIZE = max_y_id;
  GLZ_SIZE = max_z_id;
  GLYZ_SIZE = GLY_SIZE * GLZ_SIZE;
  GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;
  max_dist = 2*sqrt(GLX_SIZE*GLX_SIZE + GLY_SIZE*GLY_SIZE + GLZ_SIZE*GLZ_SIZE);

  resolution = _resolution;
  inv_resolution = 1.0 / _resolution;
  printf("map info: gl_xl %f, gl_yl %f, gl_zl %f RESOLUTION %f\n", gl_xl, gl_yl, gl_zl, resolution);
  data = new uint8_t[GLXYZ_SIZE];
  flag = new uint8_t[GLXYZ_SIZE];
  EDT = new double[GLXYZ_SIZE];
  for (int i=0;i<GLXYZ_SIZE;i++)
  {
    EDT[i] = max_dist;
  }
  memset(flag, 0, GLXYZ_SIZE * sizeof(uint8_t));
  //memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
  for (int i=0;i<GLXYZ_SIZE;i++)
  {
    data[i] = 0;
  }

  GridNodeMap = new GridNodePtr **[GLX_SIZE];
  for (int i = 0; i < GLX_SIZE; i++) {
    GridNodeMap[i] = new GridNodePtr *[GLY_SIZE];
    for (int j = 0; j < GLY_SIZE; j++) {
      GridNodeMap[i][j] = new GridNodePtr[GLZ_SIZE];
      for (int k = 0; k < GLZ_SIZE; k++) {
        Vector3i tmpIdx(i, j, k);
        Vector3d pos = gridIndex2coord(tmpIdx);
        GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
      }
    }
  }
}

void AstarPathFinder::resetGrid(GridNodePtr ptr) {
  ptr->id = 0;
  ptr->cameFrom = NULL;
  ptr->gScore = inf;
  ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids() {
  for (int i = 0; i < GLX_SIZE; i++)
    for (int j = 0; j < GLY_SIZE; j++)
      for (int k = 0; k < GLZ_SIZE; k++)
        resetGrid(GridNodeMap[i][j][k]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y,
                             const double coord_z) {
  //printf("setting obs\n");
  if (coord_x < gl_xl || coord_y < gl_yl || coord_z < gl_zl ||
      coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu)
    {/*printf("a%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n", gl_xl, gl_yl, gl_zl, gl_xu, gl_yu, gl_zu, coord_x, coord_y, coord_z);*/ return;}
  //printf("%f, %f, %f\n", coord_x, coord_y, coord_z);
  int idx_x = static_cast<int>((coord_x - gl_xl) * inv_resolution);
  int idx_y = static_cast<int>((coord_y - gl_yl) * inv_resolution);
  int idx_z = static_cast<int>((coord_z - gl_zl) * inv_resolution);
  vector<Vector3i> new_voxel;
  //printf("setting voxel: %d, %d, %d\n", idx_x, idx_y, idx_z);
  int count  = 0;
  //std::printf("%d\n", data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == uint8_t(0));
  if (idx_x == 0 || idx_y == 0 || idx_z == GLZ_SIZE || idx_x == GLX_SIZE ||
      idx_y == GLY_SIZE)
  {
    if (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == uint8_t(0)){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
      count ++;
    }
  }
  else {
    if (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == uint8_t(0)){
      printf("a\n");
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
      count ++;
    }
    if (data[(idx_x + 1) * GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] == 0){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[(idx_x + 1) * GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] = 1;
      count ++;
    }
    if (data[(idx_x + 1) * GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] == 0){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[(idx_x + 1) * GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] = 1;
      count ++;
    }
    if (data[(idx_x - 1) * GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] == 0){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[(idx_x - 1) * GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] = 1;
      count ++;
    }
    if (data[(idx_x - 1) * GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] == 0){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[(idx_x - 1) * GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] = 1;
      count ++;
    }
    if (data[(idx_x)*GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] == 0){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[(idx_x)*GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] = 1;
      count ++;
    }
    if (data[(idx_x)*GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] == 0){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[(idx_x)*GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] = 1;
      count ++;
    }
    if (data[(idx_x + 1) * GLYZ_SIZE + (idx_y)*GLZ_SIZE + idx_z] == 0){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[(idx_x + 1) * GLYZ_SIZE + (idx_y)*GLZ_SIZE + idx_z] = 1;
      count ++;
    }
    if (data[(idx_x - 1) * GLYZ_SIZE + (idx_y)*GLZ_SIZE + idx_z] == 0){
      new_voxel.push_back(Vector3i(idx_x, idx_y, idx_z));
      data[(idx_x - 1) * GLYZ_SIZE + (idx_y)*GLZ_SIZE + idx_z] = 1;
      count ++;
    }
  }
  //std::printf("%d\n", data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z]);
  updateEDT(new_voxel);
}

void AstarPathFinder::updateEDT(vector<Vector3i> new_voxel)
{
  while (!new_voxel.empty())
  {
    //std::printf("ESDF: receive %d new voxels\n", new_voxel.size());
    Vector3i pt = new_voxel.back();
    EDT[pt(0)*GLYZ_SIZE+pt(1)*GLZ_SIZE+pt(2)] = -getMinDist(pt);
    new_voxel.pop_back();
    vector<Vector3i> pt_succ;
    VoxelGetSucc(pt, pt_succ);
    
    vector<Vector3i> openSetFree;
    vector<Vector3i> openSetOccupy;
    openSetOccupy.push_back(pt);
    while (!pt_succ.empty())
    {
      Vector3i temp_pt = pt_succ.back();
      pt_succ.pop_back();
      if (isFree(temp_pt)) openSetFree.push_back(temp_pt);
    }
    while (!openSetFree.empty())
    {
      Vector3i cur_pt = openSetFree.back();
      openSetFree.pop_back();
      if (norm2(cur_pt, pt)<EDT[cur_pt(0)*GLYZ_SIZE+cur_pt(1)*GLZ_SIZE+cur_pt(2)])
      {
        //cout << "update EDT at "<< cur_pt.transpose() << "value: " << norm2(cur_pt, pt) << endl;
        EDT[cur_pt(0)*GLYZ_SIZE+cur_pt(1)*GLZ_SIZE+cur_pt(2)] = norm2(cur_pt, pt);
        //pt_succ.clear();
        VoxelGetSucc(cur_pt, pt_succ);
        while (!pt_succ.empty())
        {
          Vector3i temp_pt = pt_succ.back();
          pt_succ.pop_back();
          if (isFree(temp_pt)) openSetFree.push_back(temp_pt);
        }
      }
    }
    //std::printf("updated free space\n");
    while (!openSetOccupy.empty())
    {
      Vector3i cur_pt = openSetOccupy.back();
      openSetOccupy.pop_back();
      //pt_succ.clear();
      VoxelGetSucc(cur_pt, pt_succ);
      //printf("a\n");
      while (!pt_succ.empty())
      {
        //printf("%d\n", pt_succ.size());
        Vector3i temp_pt = pt_succ.back();
        pt_succ.pop_back();
        if (isOccupied(temp_pt))
        {
          double min_dist = getMinDist(temp_pt);
          if (min_dist > -EDT[temp_pt(0)*GLYZ_SIZE+temp_pt(1)*GLZ_SIZE+temp_pt(2)])
          {
            EDT[temp_pt(0)*GLYZ_SIZE+temp_pt(1)*GLZ_SIZE+temp_pt(2)] = -min_dist;
            //cout << "update EDT at "<< temp_pt.transpose() << "value: " << -1*min_dist << endl;
            openSetOccupy.push_back(temp_pt);
          }
        }
      }
    }
    //std::printf("updated occupied space\n");
  }

}

double AstarPathFinder::getMinDist(Eigen::Vector3i & index)
{
  if (isFree(index)) return 0.0;
  double temp_min_dist = max_dist;
  memset(flag, 0, GLXYZ_SIZE * sizeof(uint8_t));
  vector<Vector3i> occupySet;
  vector<Vector3i> succSet;
  VoxelGetSucc(index, succSet);

  while (!succSet.empty())
  {
    //printf("c\n");
    Vector3i temp_pt = succSet.back();
    succSet.pop_back();
    if (isOccupied(temp_pt))
    {
      occupySet.push_back(temp_pt);
      flag[temp_pt(0)*GLYZ_SIZE+temp_pt(1)*GLZ_SIZE+temp_pt(2)] = 1;
    }
    else
    {
      temp_min_dist = norm2(temp_pt, index);
    }
  }
  while (!occupySet.empty())
  {
    //succSet.clear();
    //printf("d\n");
    Vector3i this_pt = occupySet.back();
    occupySet.pop_back();
    if (norm2(this_pt, index)>=temp_min_dist)
      continue;
    VoxelGetSucc(this_pt, succSet);
    while(!succSet.empty())
    {
      Vector3i temp_pt = succSet.back();
      succSet.pop_back();
      if (flag[temp_pt(0)*GLYZ_SIZE+temp_pt(1)*GLZ_SIZE+temp_pt(2)]==1)
       continue;
      if (isOccupied(temp_pt))
      {
        occupySet.push_back(temp_pt);
        flag[temp_pt(0)*GLYZ_SIZE+temp_pt(1)*GLZ_SIZE+temp_pt(2)] = 1;
        //cout<<temp_pt<<endl;
      }
      else
      {
        if (norm2(temp_pt, index)<temp_min_dist)
          temp_min_dist = norm2(temp_pt, index);
      }
    }
  }
  return temp_min_dist;
}

vector<Vector3d> AstarPathFinder::getVisitedNodes() {
  vector<Vector3d> visited_nodes;
  for (int i = 0; i < GLX_SIZE; i++)
    for (int j = 0; j < GLY_SIZE; j++)
      for (int k = 0; k < GLZ_SIZE; k++) {
        // if(GridNodeMap[i][j][k]->id != 0) // visualize all nodes in open and
        // close list
        if (GridNodeMap[i][j][k]->id ==
            -1) // visualize nodes in close list only
          visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
      }

  ROS_WARN("visited_nodes size : %d", visited_nodes.size());
  return visited_nodes;
}

Vector3d AstarPathFinder::gridIndex2coord(const Vector3i &index) {
  Vector3d pt;

  pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
  pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
  pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

  return pt;
}

Vector3i AstarPathFinder::coord2gridIndex(const Vector3d &pt) {
  Vector3i idx;
  idx << min(max(int((pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
      min(max(int((pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
      min(max(int((pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);

  return idx;
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d &coord) {
  return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector3i &index) const {
  return isOccupied(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i &index) const {
  return isFree(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isOccupied(const int &idx_x, const int &idx_y,
                                        const int &idx_z) const {
  return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE &&
          idx_z >= 0 && idx_z < GLZ_SIZE &&
          (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));
}

inline bool AstarPathFinder::isFree(const int &idx_x, const int &idx_y,
                                    const int &idx_z) const {
  return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE &&
          idx_z >= 0 && idx_z < GLZ_SIZE &&
          (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr,
                                          vector<GridNodePtr> &neighborPtrSets,
                                          vector<double> &edgeCostSets) {
  neighborPtrSets.clear();
  edgeCostSets.clear();
  Vector3i neighborIdx;
  for (int dx = -1; dx < 2; dx++) {
    for (int dy = -1; dy < 2; dy++) {
      for (int dz = -1; dz < 2; dz++) {

        if (dx == 0 && dy == 0 && dz == 0)
          continue;

        neighborIdx(0) = (currentPtr->index)(0) + dx;
        neighborIdx(1) = (currentPtr->index)(1) + dy;
        neighborIdx(2) = (currentPtr->index)(2) + dz;

        if (neighborIdx(0) < 0 || neighborIdx(0) >= GLX_SIZE ||
            neighborIdx(1) < 0 || neighborIdx(1) >= GLY_SIZE ||
            neighborIdx(2) < 0 || neighborIdx(2) >= GLZ_SIZE) {
          continue;
        }

        neighborPtrSets.push_back(
            GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)]);
        edgeCostSets.push_back(sqrt(dx * dx + dy * dy + dz * dz));
      }
    }
  }
}

inline void AstarPathFinder::VoxelGetSucc(Vector3i currentPt, std::vector<Vector3i> & neighborSets)
{
  neighborSets.clear();
  for (int dx = -1; dx < 2; dx++) {
    for (int dy = -1; dy < 2; dy++) {
      for (int dz = -1; dz < 2; dz++) {

        if (dx == 0 && dy == 0 && dz == 0)
          continue;
        Vector3i neighborIdx;
        neighborIdx(0) = currentPt(0) + dx;
        neighborIdx(1) = currentPt(1) + dy;
        neighborIdx(2) = currentPt(2) + dz;

        if (neighborIdx(0) < 0 || neighborIdx(0) >= GLX_SIZE ||
            neighborIdx(1) < 0 || neighborIdx(1) >= GLY_SIZE ||
            neighborIdx(2) < 0 || neighborIdx(2) >= GLZ_SIZE) {
          continue;
        }

        neighborSets.push_back(neighborIdx);
      }
    }
  }
}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2, GridNodePtr start_node) {
  // using digonal distance and one type of tie_breaker.
  double h;
  // Du: Implement Heuristic Function for A*
  Vector3i d_node_idx;
  d_node_idx(0) = abs((node1->index)(0) - (node2->index)(0));
  d_node_idx(1) = abs((node1->index)(1) - (node2->index)(1));
  d_node_idx(2) = abs((node1->index)(2) - (node2->index)(2));
  int* distance_sort = new int[3];
  if (d_node_idx(0)<d_node_idx(1))
  {
    distance_sort[1] = d_node_idx(0);
    distance_sort[2] = d_node_idx(1);
  } else
  {
    distance_sort[1] = d_node_idx(1);
    distance_sort[2] = d_node_idx(0);
  }
  if (d_node_idx(2)<distance_sort[1])
  {
    distance_sort[0] = d_node_idx(2);
  }
  else if(d_node_idx(2)>distance_sort[2])
  {
    distance_sort[0] = distance_sort[1];
    distance_sort[1] = distance_sort[2];
    distance_sort[2] = d_node_idx(2);
  }
  else
  {
    distance_sort[0] = distance_sort[1];
    distance_sort[1] = d_node_idx(2);
  }
  //printf("distance_sort: [%d,  %d, %d]\n", distance_sort[0], distance_sort[1], distance_sort[2]);
  h = distance_sort[0]*1.718 + (distance_sort[1]-distance_sort[0])*1.414 + (distance_sort[2] - distance_sort[1]);
  Vector3i start_to_goal_idx;
  start_to_goal_idx(0) = abs((node2->index)(0) - (start_node->index)(0));
  start_to_goal_idx(1) = abs((node2->index)(1) - (start_node->index)(1));
  start_to_goal_idx(2) = abs((node2->index)(2) - (start_node->index)(2));
  double tie_breaker = sqrt(pow(abs(start_to_goal_idx(2)*d_node_idx(1)-start_to_goal_idx(1)*d_node_idx(2)), 2) + 
                            pow(abs(start_to_goal_idx(0)*d_node_idx(2)-start_to_goal_idx(2)*d_node_idx(0)), 2) +
                            pow(abs(start_to_goal_idx(1)*d_node_idx(0)-start_to_goal_idx(0)*d_node_idx(1)), 2));
  tie_breaker = 0;
  h += 0.001*tie_breaker;
  //printf("Heuristic from: [%d, %d, %d] to [%d, %d, %d]: %f tie: %f\n", node1->index(0), node1->index(1), node1->index(2), 
  //                        node2->index(0), node2->index(1), node2->index(2), h, tie_breaker);
  // Du: Implement Heuristic Function for A*
  return h;
}

void AstarPathFinder::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt) {
  ros::Time time_1 = ros::Time::now();
  std::cout<< "generate trajectory from " << start_pt.transpose() << " to " << end_pt.transpose() << std::endl;
  // index of start_point and end_point
  Vector3i start_idx = coord2gridIndex(start_pt);
  Vector3i end_idx = coord2gridIndex(end_pt);
  goalIdx = end_idx;

  // position of start_point and end_point
  start_pt = gridIndex2coord(start_idx);
  end_pt = gridIndex2coord(end_idx);

  // Initialize the pointers of struct GridNode which represent start node and
  // goal node
  GridNodePtr startPtr = new GridNode(start_idx, start_pt);
  GridNodePtr endPtr = new GridNode(end_idx, end_pt);

  // openSet is the open_list implemented through multimap in STL library
  openSet.clear();
  // currentPtr represents the node with lowest f(n) in the open_list
  GridNodePtr currentPtr = NULL;
  GridNodePtr neighborPtr = NULL;

  // put start node in open set
  startPtr->gScore = 0;
  /**
   *
   * STEP 1.1:  finish the AstarPathFinder::getHeu
   *
   * **/
  startPtr->fScore = getHeu(startPtr, endPtr, startPtr);

  startPtr->id = 1;
  startPtr->coord = start_pt;
  openSet.insert(make_pair(startPtr->fScore, startPtr));

  /**
   *
   * STEP 1.2:  some else preparatory works which should be done before while
   * loop
   *
   * **/

  double tentative_gScore;
  vector<GridNodePtr> neighborPtrSets;
  vector<double> edgeCostSets;
  printf("enter A* main loop\n");
  /**
   *
   * STEP 1.3:  finish the loop
   *
   * **/
  int count = 0;
  while (!openSet.empty()) {
    // Du: Implement A* main loop
    /*
      openSet: multimap<double, GridNodePtr>
      startPtr, endPtr, currentPtr: GridNodePtr
    */
    // get pt with max priority
    count ++;
    //printf("loop : %d, open set size: %d\n",count ,openSet.size());
    double min_fScore = openSet.begin()->first;
    currentPtr = openSet.begin()->second;
    // move current ptr from openset to closeset
    auto mapFirst = openSet.begin();
    openSet.erase(mapFirst);
    currentPtr->id = -1;
    this->AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);
    // check if pt is goal point
    if (currentPtr->coord == endPtr->coord)
    {
      printf("Exit A* main loop\n");
      return;
    }
    else
    {
      while (!neighborPtrSets.empty())
      {
        GridNodePtr succPtr = neighborPtrSets.back();
        neighborPtrSets.pop_back();
        if (isOccupied(succPtr->index)) continue;
        double succCost = edgeCostSets.back();
        edgeCostSets.pop_back();
        if (this->isOccupied(succPtr->index) || succPtr->id == -1)
        {
          continue;
        }
        else if (succPtr->id != 1)
        {
          succPtr->cameFrom = currentPtr;
          double succgScore = currentPtr->gScore + succCost;
          double succfScore = this->getHeu(succPtr, endPtr, startPtr);
          succPtr->gScore = succgScore;
          succPtr->fScore = succfScore;
          succPtr->id = 1;
          openSet.insert(make_pair(succfScore, succPtr));
        }
      } 
    }
    // Du: Implement A* main loop
  }

  // if search fails
  ros::Time time_2 = ros::Time::now();
  if ((time_2 - time_1).toSec() > 0.1)
    ROS_WARN("Time consume in Astar path finding is %f",
             (time_2 - time_1).toSec());
}

vector<Vector3d> AstarPathFinder::getPath() {
  vector<Vector3d> path;
  vector<GridNodePtr> gridPath;

  /**
   *
   * STEP 1.4:  trace back the found path
   *
   * **/
  // Du: Implement A* traceback
  GridNodePtr goalPtr = this->GridNodeMap[this->goalIdx(0)][this->goalIdx(1)][this->goalIdx(2)];
  GridNodePtr temp = goalPtr;
  path.push_back(temp->coord);
  while (temp->cameFrom != NULL)
  {
    temp = temp->cameFrom;
    path.push_back(temp->coord);
  }
  // Du: Implement A* traceback
  printf("get A* path\n");

  return path;
}

vector<Vector3d> AstarPathFinder::pathSimplify(const vector<Vector3d> &path,
                                               double path_resolution) {
  vector<Vector3d> subPath;
  /**
   *
   * STEP 2.1:  implement the RDP algorithm
   *
   * **/
  // Du: Implement RDP algorithm
  // MOOC: Motion Planning for Mobile Robotics, by Fei Gao
  double dmax = 0;
  int index = 0;
  int end = path.size();
  Vector3d st = path.front();
  Vector3d ed = path.back();
  for (int i=1;i<(end - 1);i++)
  {
    // Heron's formula
    Vector3d point = path[i];
    double a = sqrt(pow(st(0)-ed(0),2)+pow(st(1)-ed(1),2)+pow(st(2)-ed(2),2));
    double b = sqrt(pow(st(0)-point(0),2)+pow(st(1)-point(1),2)+pow(st(2)-point(2),2));
    double c = sqrt(pow(point(0)-ed(0),2)+pow(point(1)-ed(1),2)+pow(point(2)-ed(2),2));
    double p = (a+b+c)/2;
    double S = sqrt(p*(p-a)*(p-b)*(p-c));
    double d = 2*S/a;
    if (d > dmax)
    {
      dmax = d;
      index = i;
    }
  }
  vector<Vector3d> resultlist1;
  vector<Vector3d> resultlist2;
  if (dmax > path_resolution)
  {
    vector<Vector3d> path1, path2;
    for (int u=0;u <= index;u++)
    {
      path1.push_back(path[u]);
    }
    for (int k=index;k < end;k++)
    {
      path2.push_back(path[k]);
    }
    resultlist1 = pathSimplify(path1, path_resolution);
    resultlist2 = pathSimplify(path2, path_resolution);
    int a = resultlist1.size();
    int b = resultlist2.size();
    for (int j=0;j<a;j++)
    {
      subPath.push_back(resultlist1[j]);
    }
    for (int r=0;r<b;r++)
    {
      subPath.push_back(resultlist2[r]);
    }
  } else
  {
    subPath.push_back(path.front());
    subPath.push_back(path.back());
  }
  // Du: Implement RDP algorithm
  return subPath;
}

Vector3d AstarPathFinder::getPosPoly(MatrixXd polyCoeff, int k, double t) {
  Vector3d ret;
  int _poly_num1D = (int)polyCoeff.cols() / 3;
  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * _poly_num1D, _poly_num1D);
    VectorXd time = VectorXd::Zero(_poly_num1D);

    for (int j = 0; j < _poly_num1D; j++)
      if (j == 0)
        time(j) = 1.0;
      else
        time(j) = pow(t, j);

    ret(dim) = coeff.dot(time);
    // cout << "dim:" << dim << " coeff:" << coeff << endl;
  }

  return ret;
}

int AstarPathFinder::safeCheck(MatrixXd polyCoeff, VectorXd time) {
  int unsafe_segment = -1; //-1 -> the whole trajectory is safe
  /**
   *
   * STEP 3.3:  finish the sareCheck()
   *
   * **/
  // Du: safety check of polynomial traj
  // MOOC: Motion Planning for Mobile Robots, by Fei Gao
  
  // Du: safety check of polynomial traj
  return unsafe_segment;
}