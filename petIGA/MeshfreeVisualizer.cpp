// Petsc/PetIGA header files
#include "petiga.h"
#include "../src/petigagrid.h"
#include <petscsys.h>
#include <../../../../../src/sys/fileio/mprint.h>
#include <iomanip>
#include <petscblaslapack.h>
#include <petsc/private/tsimpl.h>


// Boost parallelization libraries
#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/distributed/distributed_graph_utility.hpp>
#include <boost/graph/distributed/local_subgraph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/version.hpp>


// Standard Library includes
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
using namespace Teuchos;


#define SQ(A) ((A)*(A))
#define for1(i,n) for(int i=0; i<(n); i++)
#define for2(i,j,n) for(int i=0; i<(n); i++) for(int j=0; j<(n); j++)
#define for3(i,j,k,n) for(int i=0; i<(n); i++) for(int j=0; j<(n); j++) for(int k=0; k<(n); k++)
#define for4(i,j,k,l,n) for(int i=0; i<(n); i++) for(int j=0; j<(n); j++) for(int k=0; k<(n); k++) for(int l=0; l<(n); l++)
#define for5(i,j,k,l,m,n) for(int i=0; i<(n); i++) for(int j=0; j<(n); j++) for(int k=0; k<(n); k++) for(int l=0; l<(n); l++) for(int m=0; m<(n); m++)
#define for6(i,j,k,l,m,o,n) for(int i=0; i<(n); i++) for(int j=0; j<(n); j++) for(int k=0; k<(n); k++) for(int l=0; l<(n); l++) for(int m=0; m<(n); m++) for(int o=0; o<(n); o++)
#define for7(i,j,k,l,m,o,p,n) for(int i=0; i<(n); i++) for(int j=0; j<(n); j++) for(int k=0; k<(n); k++) for(int l=0; l<(n); l++) for(int m=0; m<(n); m++) for(int o=0; o<(n); o++) for(int p=0; p<(n); p++)
#define for8(i,j,k,l,m,o,p,q,n) for(int i=0; i<(n); i++) for(int j=0; j<(n); j++) for(int k=0; k<(n); k++) for(int l=0; l<(n); l++) for(int m=0; m<(n); m++) for(int o=0; o<(n); o++) for(int p=0; p<(n); p++) for(int q=0; q<(n); q++)

// Quantities that will be reduced into through MPI and Broadcast
int mpiSize;
double totalInitialExplosiveVolume;
double totalCurrentExplosiveVolume;
double totalExplosiveMass;
double totalExplosiveVolume;
int totalNumNodes;
int num_PD_nodes;

//// Data Structures for Immersed-IGA-PD FSI ////

//Information we need for the background discretization
//For example Lx,Ly,Lz are the lengths of the background domain in x,y,z
//Nx,Ny,Nz are the numbers of elements in each directions
//A1,Ap,An are the n+1,alpha,n levels of the acceleration
//V1,Vp,Vn are the n+1,alpha,n levels of the velocity
//dA is the increment of the acceleration. This is what we are solving for
//Alpha_m,Alpha_f,Gamma,Beta are parameters for the generalized alpha method
typedef struct {

  IGA iga;
  IGA iga_energy;
  IGA iga_strain;
  PetscReal Lx,Ly,Lz,mu,lamda,kappa,temp0,p0,R,Cp,Cv,spacing;
  Vec A1,Ap,An,V1,Vp,Vn,D1,Dp,Dn,Aa,Va,Da,V0,A0,D0,dA;
  PetscInt  Nx,Ny,Nz,max_its;
  Mat Ainv;
  Mat MFinv;

  PetscInt *processor_numElX;
  PetscInt *processor_numElY;
  PetscInt *processor_numElZ;
  PetscReal *xhat;
  PetscReal *yhat;
  PetscReal *zhat;

  PetscReal totalInitialExplosiveVolume;
  PetscReal totalCurrentExplosiveVolume;
  PetscReal totalExplosiveMass;
  PetscReal H[3];
  PetscReal TimeRestart;
  PetscInt  numFluidNodes;

  PetscInt  StepRestart;
  PetscInt  stepNumber;
  PetscInt  it;
  PetscInt  FreqRestarts;
  PetscInt xdivint;
  PetscInt ydivint;
  PetscInt zdivint;
  PetscInt nen;
  PetscReal *GlobalForces;
  PetscReal *GlobalVel;
  PetscReal *GlobalDisp;
  PetscReal *COORD;
  //PD Nodes on which the boundary conditions are enforeced. These will be read from node_list files
  //ID's in these files correspond to fd PD_ID's, so anytime a bc needs to be set, the fd.PD_ID and the PD_ID_BC
  //can be checked.
  PetscInt  *PD_ID_BC;

  PetscReal thickness;
  PetscReal Alpha_m,Alpha_f,Gamma,Beta;
  PetscInt num_Owned_Points;

} AppCtx;

// info to identify a vertex; need to be able to get unique ID on the fly
// w/out communication
class VertexID{
public:
  // three-part ID, using current processor rank, the rank of
  // the processor that created the vertex, and the local ID.  birthCert
  // and localID uniquely identify the vertex.  includeing current rank
  // in the VertexID and then returning it as a hash is the only way
  // I could find to explicitly control which process the vertex goes to.
  int rank;
  int birthCert;
  int localID;

  // in calls to constructors, birthCert_ should ALWAYS be the calling
  // task's MPI rank.  rank_ can be any other task, and determines
  // what task will own the new vertex created with this VertexID

  // constructor: updates local count automatically
  VertexID(int rank_, int birthCert_, int *localCount){
    rank = rank_;
    birthCert = birthCert_;
    localID = (*localCount);
    (*localCount)++;
  } // end constructor

  // constructor: does not update local count
  VertexID(int rank_, int birthCert_, int local_){
    rank = rank_;
    birthCert = birthCert_;
    localID = local_;
  } // end non-updating constructor

  // default constructor
  VertexID(){}

  void doCopy(const VertexID &id){
    rank = id.rank;
    birthCert = id.birthCert;
    localID = id.localID;
  } // end doCopy

  // destructor
  ~VertexID(){}

  // assignment
  VertexID &operator=(const VertexID &id){
    doCopy(id);
  } // end assignment

  template<typename Archiver>
  void serialize(Archiver& ar, const unsigned int version){
    ar & rank & birthCert & localID;
  } // end serialize
}; // end VertexID

BOOST_IS_MPI_DATATYPE(VertexID);
// Particle information (this could be part of FieldData instead, but the idea
// is that this is "small" and contains less stuff that everything we
// need for computations at a particle, e.g., things used in decisions
// about graph connectivity, iteration, activity/inactivity, etc.)
class ParticleInfo{
public:
  // REMEMBER: if new field added, also add to serialize method
  double currentCoord[3];
  double initialCoord[3];
  double tempCoord[3];
  bool   isTask;

  // TODO: other flags indicating stuff about the vertex

  ParticleInfo(){
    for1(i,3)
      currentCoord[i] = 0.0;
    isTask = false;
  } // end construtor

  ~ParticleInfo(){}

  template<typename Archiver>
  void serialize(Archiver& ar, const unsigned int version){
    ar & currentCoord & initialCoord & tempCoord & isTask;
  } // end serialize
}; // end class ParticleInfo

BOOST_IS_MPI_DATATYPE(ParticleInfo);

// Holds state information for a particle; idea is that this is larger in
// memory than ParticleInfo, and holds all sorts of tensors, etc., that are
// only really needed in constitutive routines
class FieldData{
public:
  // REMEMBER: if new field added, also add to serialize method

  // Variables associated with assembly of foreground PD residual and
  // Peridigm integration
  double inertia[3];
  double residual[3];
  double internalForce[3];
  double bodyForce[3];
  double referenceNodalVolume;
  double alphaNodalVolume;
  double alphaNodalDensity;
  double referenceDensity;
  //
  int material;

  double nodalVolume;

  //For quadrature volume update, new variables:
  double nodalDensity;
  double nodalDensityInitial;
  double nodalPressure;
  double nodalVolumeInitial;

  int    Boundary;
  double totalPhysicalDisplacement[3];
  double totalPhysicalVelocity[3];
  double totalStrain[6];
  double totalStress[6];
  double totalStress0[6];
  double totalStrain0[6];
  double ductile_threshold0;
  double ductile_threshold;
  double brittle_threshold0;
  double brittle_threshold;
  double damage;
  double damage0;
  double currentDeformationGradient[9];
  double DeformationGradientOld[9];
  double alphaDeformationGradient[9];
  double velocityGradient[9];
  double determinantCurrentDeformationGradient;
  double determinantAlphaDeformationGradient;

  double totalPhysicalDisplacementOldIteration[3];
  double totalPhysicalDisplacementOldStep[3];
  double totalPhysicalVelocityOldIteration[3];
  double totalPhysicalVelocityOldStep[3];
  double totalPhysicalAcceleration[3];
  double totalPhysicalAccelerationOldIteration[3];
  double totalPhysicalAccelerationOldStep[3];
  double AccelerationIncrement[3];
  double currentCoord[3];
  double referenceCoord[3];


  double ductile_energy;
  double brittle_energy;
  int    Inside; //Flag to state if the particle is inside the background computational domain. 1 for inside, 0 for outside
  int    ID;
  int    ID_PD;

  double effectiveStrainRate;

  int flag;
  int flag0;

  /// ### Penalty Associated Variables ### ///
  double interpolatedVelocity[3];
  double penaltyParameter;
  double referencePenaltyParameterInternal;
  double referencePenaltyParameterInertia;
  bool   flyingPoint;
  double penaltyForce[3];
  //////////////////////////////////////////

  FieldData(){}
  ~FieldData(){}

  template<typename Archiver>
  void serialize(Archiver& ar, const unsigned int version){
    ar & interpolatedVelocity & penaltyParameter & referencePenaltyParameterInternal & referencePenaltyParameterInertia & flyingPoint & penaltyForce & nodalVolume & nodalDensity & nodalDensityInitial & nodalPressure & nodalVolumeInitial & Boundary & totalPhysicalDisplacement & totalPhysicalVelocity & totalStrain & totalStress & totalStrain0 & totalStress0
  & ductile_threshold0 & ductile_threshold & brittle_threshold0 & brittle_threshold & damage0 & damage & currentDeformationGradient & velocityGradient
  & determinantCurrentDeformationGradient & totalPhysicalDisplacementOldIteration & totalPhysicalDisplacementOldStep & totalPhysicalVelocityOldIteration
  & totalPhysicalVelocityOldStep & totalPhysicalAcceleration & totalPhysicalAccelerationOldIteration & totalPhysicalAccelerationOldStep & AccelerationIncrement
    & ductile_energy & brittle_energy & Inside & DeformationGradientOld & alphaDeformationGradient & determinantAlphaDeformationGradient & effectiveStrainRate
  & flag & flag0 & inertia & residual & bodyForce & internalForce & referenceNodalVolume & alphaNodalDensity & alphaNodalVolume & ID & ID_PD & material & referenceDensity & referenceCoord & currentCoord;
  } // end serialize
}; // end class fieldData

BOOST_IS_MPI_DATATYPE(FieldData);

// class to hold all of the stuff at a graph vertex
class VertexData{
public:

  // the approach i've taken is to store parts of data in other classes, then
  // have instances of those classes as members.  the reason i do this instead
  // of storing everything as members of VertexData is that, to access
  // members of VertexData from vertex descriptors (typedef'd as "Vertex")
  // we'd need a new property map for each one.  this would make
  // adding/removing new fields/data associated with particles very tedious.

  // This exists to uniquely-identify a vertex
  VertexID id;

  // this stores a small amount of info used by things other than
  // constitutive routines
  ParticleInfo info;

  // this stores all of the physical data associated with a point (except
  // for its location, which is used in determining graph connectivity,
  // and is stored in info)
  FieldData fd;

  // given ID only
  VertexData(VertexID id_){
    id = id_;
    info = ParticleInfo();
    fd = FieldData();
  } // end ID only constructor

  // given ID, field data
  VertexData(VertexID id_, ParticleInfo info_){
    id = id_;
    info = info_;
    fd = FieldData();
  } // end ID only constructor

  VertexData(){}
  ~VertexData(){}

  template<typename Archiver>
  void serialize(Archiver& ar, const unsigned int version){
    ar & id & fd & info;
  } // end serialize

}; // end VertexData

BOOST_IS_MPI_DATATYPE(VertexData);

// hash function for VertexID class
// appaerently this determines which process a named vertex goes to
// (you'd think the docs would say this somewhere...)
std::size_t hash_value(const VertexID &v) {
  //return boost::hash_value((size_t)v.localID + (size_t)v.rank);
  return (size_t)v.rank + ((size_t)mpiSize)*((size_t)v.localID);
  //return (size_t)v.rank;
} // end hash function

// overload operator outside of class?
bool operator==(const VertexID &id1, const VertexID &id2){
  return (id1.rank==id2.rank)&&(id1.localID==id2.localID)
    &&(id1.birthCert==id2.birthCert);
} // end overloaded ==

// This snippet is used to allow acces of vertices through their VertexID
// attributes
namespace boost { namespace graph {
    template<>
    struct internal_vertex_name<VertexData>
    {
      typedef multi_index::member<VertexData, VertexID, &VertexData::id> type;
    };
  }
} // end activating named vertices

// This specifies what to do if a vertex is requested via ID, but no such
// vertex exists
namespace boost { namespace graph {
    template<>
    struct internal_vertex_constructor<VertexData>
    {
      typedef vertex_from_name<VertexData> type;
    };
  }
} // end allowing new verticies

using namespace boost;
using boost::graph::distributed::mpi_process_group;
// there is an apparent bug in the parallel boost graph library that causes
// calls to the existing remove_vertex() function to cause
// compiler errors, so I introduced my own version that comments out the
// offending line.  No idea whether this could potentially make trouble, but
// it seems to work okay...
template<PBGL_DISTRIB_ADJLIST_TEMPLATE_PARMS>
void
remove_vertex_dk(typename PBGL_DISTRIB_ADJLIST_TYPE::vertex_descriptor u,
        PBGL_DISTRIB_ADJLIST_TYPE& graph)
{
  typedef typename PBGL_DISTRIB_ADJLIST_TYPE::graph_type graph_type;
  typedef typename graph_type::named_graph_mixin named_graph_mixin;
  BOOST_ASSERT(u.owner == graph.processor());
  static_cast<named_graph_mixin&>(static_cast<graph_type&>(graph))
    .removing_vertex(u, boost::graph_detail::iterator_stability
         (graph.base().m_vertices));

  // this line caused an error ///////
  //g.distribution().clear();
  ////////////////////////////////////

  remove_vertex(u.local, graph.base());
}

// bidirectional graph, using a list for the underlying data structure
typedef adjacency_list<listS, // out edge list type
           distributedS<mpi_process_group, listS>, // vertex list
           //bidirectionalS,
           directedS, // graph type
           VertexData, // vertex property
           no_property, // edge property
           no_property, // graph property
           listS> // edge list
Graph;

// renaming convoluted boost types for easier use
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::vertex_iterator VertexIterator;
typedef graph_traits<Graph>::adjacency_iterator AdjIterator;
typedef graph_traits<Graph>::edge_iterator EdgeIterator;
typedef graph_traits<Graph>::out_edge_iterator OutEdgeIterator;

// these allow us to access members of the vertex property type VertexData
// given a vertex descriptor
property_map<Graph, VertexID VertexData::*>::type
id_property;
property_map<Graph, ParticleInfo VertexData::*>::type
info_property;
property_map<Graph, FieldData VertexData::*>::type
fd_property;

// class that wraps around Graph and provides access to property maps
// and methods for refinement and coarsening
class ParticleManager{
public:
  Graph graph;
  int localVertexCounter;
  int myRank;
  int nProc;
  Vertex myTaskVertex;
  vector<Vertex> taskVertexDescriptors;

  // copy constructor
  ParticleManager(const ParticleManager &manager){}

  void generatePropertyMaps(){
    id_property = get(&VertexData::id, graph);
    id_property.set_consistency_model
      (boost::parallel::consistency_model::cm_bidirectional);
    info_property = get(&VertexData::info, graph);
    info_property.set_consistency_model
      (boost::parallel::consistency_model::cm_bidirectional);
    fd_property = get(&VertexData::fd, graph);
    fd_property.set_consistency_model
      (boost::parallel::consistency_model::cm_bidirectional);
  } // end generatePropertyMaps

  // synchronize processors
  void sync(){
    synchronize(graph.process_group());
  } // end sync

  // constructor
  ParticleManager(){

    generatePropertyMaps();
    myRank = graph.process_group().rank;
    nProc = graph.process_group().size;
    taskVertexDescriptors = vector<Vertex>();

    // add graph vertices corresponding to processor subdomains
    for1(rank,nProc){
      // dummy local vertex counter
      int dummyCounter = 0;
      // create VertexID that hashes to rank and has local vert number or zero
      VertexID id = VertexID(rank,rank,&(dummyCounter));
      // create vertex from VertexID
      VertexData vd = VertexData(id);
      vd.info.isTask = true;
      // actually add vertex to graph if mine
      if(rank == myRank)
  // keep track of my task vertex descriptor
  myTaskVertex = add_vertex(vd,graph);
    } // rank
    sync();

    // add descriptors to task vertexes
    for1(rank,nProc){
      // dummy local vertex counter
      int dummyCounter = 0;
      // create VertexID that hashes to rank and has local vert number or zero
      VertexID id = VertexID(rank,rank,&(dummyCounter));
      // add corresponding descriptor to list of task verts on all ranks
      taskVertexDescriptors.push_back(*find_vertex(id,graph));
    } // rank
    sync();

    // set local counter to 1 to reflect task vert added
    localVertexCounter = 1;

  } // end constructor


  // obtain, in O(1) time, the rank of the processor whose subdomain contains
  // the point
  int pointRank(double *x, AppCtx *user){

    PetscInt elemX,elemY,elemZ;
    int i;

    //This has to be redone for non-uniform meshes:
    elemX = (int) (x[0]/user->H[0]);
    elemY = (int) (x[1]/user->H[1]);
    elemZ = (int) (x[2]/user->H[2]);
    if (elemX == user->iga->elem_sizes[0] && x[0] <= user->Lx){
      elemX = user->iga->elem_sizes[0]-1;
    }
    if(elemY == user->iga->elem_sizes[1] && x[1] <= user->Ly){
      elemY = user->iga->elem_sizes[1]-1;
    }
    if(elemZ == user->iga->elem_sizes[2] && x[2] <= user->Lz){
      elemZ = user->iga->elem_sizes[2]-1;
    }


    PetscInt rankX=-1,rankY=-1,rankZ=-1;
    PetscInt check = 0;
    for (i=0;i<user->iga->proc_sizes[0];i++){
      check += user->processor_numElX[i];
      if (elemX<check){
        rankX = i;
        break;
      }
    }
    check = 0;
    for (i=0;i<user->iga->proc_sizes[1];i++){
      check += user->processor_numElY[i];
      if (elemY<check){
        rankY = i;
        break;
      }
    }
    check = 0;
    for (i=0;i<user->iga->proc_sizes[2];i++){
      check += user->processor_numElZ[i];
      if (elemZ<check){
        rankZ = i;
        break;
      }
    }


    PetscInt myRank = rankZ*user->iga->proc_sizes[0]*user->iga->proc_sizes[1] + rankY*user->iga->proc_sizes[0]+rankX;
    if ((elemX < 0) || (elemY < 0) || (elemZ < 0)){
      myRank = -1;
      PetscPrintf(PETSC_COMM_SELF, "Error! Point at x = %e %e %e returned negative Element values: elem = %d %d %d with H = %e %e %e\n", x[0], x[1], x[2], elemX, elemY, elemZ, user->H[0], user->H[1], user->H[2]);
      exit(0);
    }
    if ((elemX >= user->iga->elem_sizes[0]) || (elemY >= user->iga->elem_sizes[1]) || (elemZ >= user->iga->elem_sizes[2])){
      PetscPrintf(PETSC_COMM_SELF, "Error! Point at x = %e %e %e returned Element values which exceed the number of elements: elem = %d %d %d with H = %e %e %e\n", x[0], x[1], x[2], elemX, elemY, elemZ, user->H[0], user->H[1], user->H[2]);
      myRank = -1;
      exit(0);
    }
    return myRank;
  } // end pointRank

  // connect particles to their corresponding task vertexes
  // "reconnect" is whether or not verts have previously been connected
  void connectVertsToTasks(bool reconnect,AppCtx *user){

    // vector to store information about what edges to add
    vector<pair<Vertex,Vertex>> edgesToAdd = vector<pair<Vertex,Vertex>>();

    // if we're just updating existing connectivity:
    if(reconnect){
      // iterate over all out-edges of this task's vertex.  delete ones
      // that are no longer valid, and make a note (in edgesToAdd) to add
      // any new edges that are needed.
 pair<OutEdgeIterator,OutEdgeIterator> its
  = out_edges(myTaskVertex,graph);
      auto it=its.first;
      while(it != its.second){
  Edge edge = *it;
  // vertex descriptor for particle on the other end of the edge
  Vertex v = target(edge,graph);
  // get the vertex's position, and use it to decide whether or not it
  // is still located in this processor's subdomain
  ParticleInfo info = get(info_property,v);
  int pr = pointRank(&info.currentCoord[0],user);
  // if the particle is no longer located on this processor's subdomain
  if(pr != myRank){
    // make a note to add a new edge from the task whose subdomain
    // the particle is located in.
          if (pr >= 0){
    edgesToAdd.push_back
      (pair<Vertex,Vertex>(taskVertexDescriptors[pr],v));
    }
    // remove the old edge:
    // note: removing an edge will invalidate the iterators pointing to
    // it, so we need to make a copy, advance the original, then
    // pass the copy to remove_edge()
    auto itTemp = it;
    ++it;
    // the version taking an iterator (instead of *itTemp) is faster
    // for directed graphs
    remove_edge(itTemp,graph);
  }else{
    // if the particle is still in this processor's subdomain, just
    // move on to the next edge
    ++it;
  } // end if
      } // it
    }else{ // otherwise, assume we're generating connectivity from scratch

      // remove existing edges (assume O(E/V) time...?)
      // if we're really only calling this on the first step, we don't need
      // to call clear_...(), but this branch is useful for debugging.
      clear_out_edges(myTaskVertex,graph);

      // loop over all vertices in the graph, determine which processor's
      // subdomain they are in, then make a note to add the corresponding
      // edge.
      BGL_FORALL_VERTICES(v,graph,Graph){
         ParticleInfo info = get(info_property,v);
           if(!info.isTask){
    // add edge from task to vertex
    edgesToAdd.push_back
      (pair<Vertex,Vertex>
       (taskVertexDescriptors[pointRank(&info.currentCoord[0],user)],v));
  } // end if is a particle, not a task
      } // v
    } // end switch on reconnecting
    sync();

    // actually create new edges from information stored in edgesToAdd

    for1(i,edgesToAdd.size()){
      add_edge(edgesToAdd[i].first,edgesToAdd[i].second,graph);
    } // i
    sync();

    // generatePropertyMaps();
    // sync();
    //PetscPrintf(PETSC_COMM_WORLD,"Return connectVertsToTasks");
  } // end connectVertsToTasks
}; // end ParticleManager class

#define SQ(A) ((A)*(A))
typedef struct
{
    PetscReal initialTime;
    PetscReal finalTime;
    PetscReal currentTime;
    PetscReal gamma;
    PetscReal beta;
    PetscReal timeStep;
    PetscInt  stepNumber;

    PetscReal max_compressive_strain;
    PetscReal max_tensile_strain;
    PetscReal max_compressive_stress;
    PetscReal max_tensile_stress;
    PetscReal brittle_initial_threshold;
    PetscReal ductile_initial_threshold;

    PetscReal   youngModulus;
    PetscReal   poissonRatio;
    PetscReal   density;
    PetscReal   lambda;
    PetscReal   mu;
    PetscReal   massDamping;

    PetscReal   SigmaYinitial;  //Plasticity

    PetscInt    FreqResults;
    PetscInt    numPoints;
    PetscInt    numNodes;
    PetscReal   Alpha_f;

    PetscReal densityRDX;

    PetscInt rateEffects;

    PetscReal penaltyConstant;
    PetscBool DamageModeling;
    PetscReal damageCriticalStress;
    PetscReal damageCriticalEpsilonPlastic;
    PetscReal thresholdDamageForPenalty;


} PARAMETERS;

typedef struct {
  PetscScalar rho,ux,uy,uz,temp;
} Field;
//// End Data structures for immersed-IGA ////

//// Re-definition of PetIGA functions ////
#undef  __FUNCT__
#define __FUNCT__ "IGALocateElement_1"
PetscBool IGALocateElement_1(IGA iga,PetscReal *pnt,IGAElement element)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscInt i,j,e,m,dim=iga->dim;
  PetscInt *ID = element->ID;
  PetscInt *width = element->width;
  PetscInt *start = element->start;
  PetscScalar *U;

  element->nen  = 1;
  element->nval = 0;
  element->nvec = 0;
  element->nmat = 0;

  for(i=0;i<dim;i++){
    element->nen *= (iga->axis[i]->p+1);
    U = iga->axis[i]->U;
    ID[i] = 0;
    PetscReal h = U[iga->axis[i]->p+1] - U[iga->axis[i]->p];
    PetscReal deltau = pnt[i] - U[0];
    e = (PetscInt) (deltau/h);
    ID[i] = e;
    /* find which nonzero span this point is located in */
   // for(j=0;j<m;j++){
      //if(U[j+1]-U[j]>1.0e-13) e += 1;
      //if(pnt[i] > U[j] && pnt[i] <= U[j+1]) ID[i] = e;
    //}
    /* reject if the element is not in this partition */
     if(ID[i] < iga->elem_start[i] || ID[i] >= iga->elem_start[i]+iga->elem_width[i]) return PETSC_FALSE;
  }

  element->index = 0;
{
    PetscErrorCode ierr;
    ierr = IGAElementBuildClosure(element);CHKERRCONTINUE(ierr);
    if (PetscUnlikely(ierr)) PetscFunctionReturn(PETSC_FALSE);
    ierr = IGAElementBuildFix(element);CHKERRCONTINUE(ierr);
    if (PetscUnlikely(ierr)) PetscFunctionReturn(PETSC_FALSE);
  }

  return PETSC_TRUE;
}




PETSC_STATIC_INLINE
PetscBool IGAElementNextFormIJacobian(IGAElement element,IGAFormIJacobian *jac,void **ctx)
{
  IGAForm form = element->parent->form;
  if (!IGAElementNextForm(element,form->visit)) return PETSC_FALSE;
  *jac = form->ops->IJacobian;
  *ctx = form->ops->IJacCtx;
  return PETSC_TRUE;
}


//// End Re-Definition of PetIGA functions ////


//// Template functions ////

#undef __FUNCT__
#define __FUNCT__ "FG_Template"
PetscErrorCode FG_Template(AppCtx *user, PARAMETERS *par, ParticleManager &manager)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;

  PetscReal pt[3];
  pair<OutEdgeIterator,OutEdgeIterator> its = out_edges(manager.myTaskVertex,manager.graph);

  IGAProbe p;
  IGAProbeCreate(user->iga,user->V1,&p);
  IGAProbeSetCollective(p, PETSC_FALSE);

  for(auto it=its.first; it != its.second; ++it){
    Edge edge = *it;
    Vertex v = target(edge,manager.graph);
    ParticleInfo info = get(info_property,v);
    FieldData fd = get(fd_property,v);
    pt[0] =  info.currentCoord[0]/user->Lx;
    pt[1] =  info.currentCoord[1]/user->Ly;
    pt[2] =  info.currentCoord[2]/user->Lz;
    ierr = IGAProbeSetPoint(p,pt);CHKERRQ(ierr);
    put(fd_property,v,fd);
  }

PetscFunctionReturn(0);
}

////////////////////////////

//// I/O ////
#undef  __FUNCT__
#define __FUNCT__ "input"
// Reads foreground input files with format Pre-Processor V3:
//This function reads the input files for the particle data, specifically foreground(pr+1).dat
// Read file structured as:
// #nodes 0 0 0 0 0 0 0 0
// ID x y z dx dy dz material_ID vol
PetscErrorCode input (PARAMETERS *par,ParticleManager &manager, AppCtx *user)
{

  PetscFunctionBegin;
  //PetscPrintf(PETSC_COMM_WORLD," Input Begin");
  PetscErrorCode  ierr;
  PetscInt i,j;
  PetscInt num_proc;
  PetscInt temp;
  PetscReal tempR;
  int pr;

  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //PetscPrintf(PETSC_COMM_WORLD," Processor number: %d", rank);


  ostringstream convert;
  convert << "foreground" << rank+1 << ".dat";
  string fname = convert.str();
  ifstream fin;
  fin.open(fname.c_str());

  fin >> par->numNodes;
  for (j = 0; j<6; j++){
    fin >> tempR; //Skip over the padding zeros in first line (8 zeros)
  }
  fin >> temp; // Material_ID field is an integer
  fin >> tempR;


  par->numPoints = par->numNodes;
 //PetscPrintf(PETSC_COMM_SELF," %d \n",par->numNodes);
  // add a bunch of verticies
            for1(i,par->numNodes){
              ParticleInfo info = ParticleInfo();
              FieldData fd = FieldData();
              fin >> fd.ID;//ID
              fin >> info.currentCoord[0];//x
              fin >> info.currentCoord[1];//y
              fin >> info.currentCoord[2];//z
              fin >> tempR;//x-width
              fin >> tempR;//y-width
              fin >> tempR;//z-width

              //Material Flag:
              fin >> fd.material;
              // Current definition:
              // 0 - Composite Shell
              // 1 - RDX
              // 2 - Air

              //initialize NodalVolume = x*y*z
              fin >> fd.nodalVolume;

              //Initialize reference and current quantities:
              info.initialCoord[0] = info.currentCoord[0];
              info.initialCoord[1] = info.currentCoord[1];
              info.initialCoord[2] = info.currentCoord[2];
              fd.nodalVolumeInitial = 0.0;
              fd.nodalVolumeInitial += fd.nodalVolume;

              //Debugging
              //PetscPrintf(PETSC_COMM_WORLD,"volI %e vol %e \n",fd.nodalVolume, fd.nodalVolumeInitial);


              fd.Boundary = 0;
              fd.Inside = 1;

              //Find PD solid boundaries
              if(fd.material==0){
                fd.damage = 0;
                fd.flyingPoint = false;
                // if(info.initialCoord[0] < 1e-15+user->spacing/2.0+0.051/2.0 || info.initialCoord[1] < 1e-15+user->spacing/2.0+0.051/2.0 || abs(info.initialCoord[0]-user->Lx+user->spacing/2.0)<1e-15+0.051/2.0 || abs(info.initialCoord[1]+user->spacing/2.0-(user->Ly))<1e-15+0.051/2.0){
                //   fd.Boundary=1;
                // }
                // if(info.initialCoord[0] < 1e-15+user->spacing/2.0 || info.initialCoord[1] < 1e-15+user->spacing/2.0 || abs(info.initialCoord[0]-user->Lx+user->spacing/2.0)<1e-15 || abs(info.initialCoord[1]+user->spacing/2.0-(user->Ly))<1e-15){
                //   fd.Boundary=1;
                // }
                //fd.nodalVolume = fd.nodalVolume*user->thickness;
                PetscReal meshSize = 0.305/150.0;//(par->puntos[i].support[0] + par->puntos[i].support[1])/(2.0*par->supportFactor); // average particle spacing
                fd.referencePenaltyParameterInternal = par->penaltyConstant * par->youngModulus * par->timeStep / (meshSize * meshSize);
                fd.referencePenaltyParameterInertia = par->penaltyConstant * par->density / par->timeStep;
                fd.penaltyParameter = fd.referencePenaltyParameterInternal;
              }

              for (j=0;j<6; j++){
                fd.totalStrain[j]  = 0.0;
                fd.totalStrain0[j] = 0.0;
                fd.totalStress[j]  = 0.0;
                fd.totalStress0[j] = 0.0;
              }
              for (j=0;j<3; j++){
                fd.totalPhysicalDisplacement[j] = 0.0;
                fd.totalPhysicalVelocity[j]     = 0.0;
                info.currentCoord[0] = info.currentCoord[0]+0.00000;
                info.currentCoord[1] = info.currentCoord[1]+0.00000;
                info.currentCoord[2] = info.currentCoord[2]+0.00000;
              }
              if (par->stepNumber == 0){
                    pr = manager.pointRank(&info.currentCoord[0],user);
                    VertexData vd = VertexData(VertexID(pr,rank,&manager.localVertexCounter),info);
                    vd.fd = fd;
                    add_vertex(vd,manager.graph);
              }
            }

          	user->nen  = 0;
          	user->nen  = user->iga->axis[0]->p + 1;
          	user->nen *= user->iga->axis[1]->p + 1;
          	user->nen *= user->iga->axis[2]->p + 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormInitialCondition"
PetscErrorCode FormInitialCondition(IGA iga,PetscReal t,Vec U,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  DM da;
  PetscInt dof = iga->dof;
  PetscInt dim = iga->dim;

  ierr = IGACreateNodeDM(iga,dof,&da);CHKERRQ(ierr);
  Field ***u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;

  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

  PetscInt i,j,k;
  PetscInt nodesX  = iga->geom_lwidth[0], nodesY  = iga->geom_lwidth[1], nodesZ  = iga->geom_lwidth[2];
  PetscInt gnodesX = iga->geom_gwidth[0], gnodesY = iga->geom_gwidth[1];
  PetscReal hx = 0.0;

  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      for(k=info.zs;k<info.zs+info.zm;k++){
      u[k][j][i].ux   =  0.0;
      u[k][j][i].uz   =  0.0;
      u[k][j][i].uy   =  0.0;
      u[k][j][i].rho  =  1.0;
      u[k][j][i].temp =  1.0;
    }
  }
}
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "ParticleDistribute"
PetscErrorCode ParticleDistribute(PARAMETERS *par, AppCtx *user, ParticleManager &manager)
{
  // Based on the information about the PD initialization, update ID, ID_PD and material info
  // so that other quantities in the iteration loop can be computed conditioned on this information.
  //Debugging
  //PetscPrintf(PETSC_COMM_SELF, "%d\n", user->numFluidNodes);

  PetscFunctionBegin;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  user->totalInitialExplosiveVolume = 0.0;
  user->totalCurrentExplosiveVolume = 0.0;
  pair<OutEdgeIterator,OutEdgeIterator> its = out_edges(manager.myTaskVertex,manager.graph);
  // Here we map to the PD object ID so that the kinematic fields can be updated
  for(auto it=its.first; it != its.second; ++it){
    Edge edge = *it;
    Vertex v = target(edge,manager.graph);
    FieldData fd = get(fd_property,v);
    ParticleInfo info = get(info_property,v);
    fd.ID_PD = -1;
  if(fd.material==0 && fd.ID>user->numFluidNodes && !info.isTask){
    fd.ID_PD = fd.ID-user->numFluidNodes-1;
    //Debugging
    PetscReal x = info.currentCoord[0];
    PetscReal y = info.currentCoord[1];
    PetscReal z = info.currentCoord[2];

    if(fd.ID_PD<0){
    PetscPrintf(PETSC_COMM_SELF, "PD_ID < 0 : number of processors < number of input files \n");
    exit(0);
    }
  }

  // In Pre-Processor, index starts at 1, so we correct so that it is 0->N-1 so it can be used for
  // indexing arrays which correspond to i = 1 : N
  if(fd.material==1){
  user->totalInitialExplosiveVolume+=fd.nodalVolume;
  user->totalCurrentExplosiveVolume+=fd.nodalVolume;
  }
  put(fd_property,v,fd);

}

  //Debugging
  // if(rank==0 || rank==1){
  //   PetscPrintf(PETSC_COMM_SELF, "Number of air particles on rank %d = %d and PD Particles = %d and nodes = %d\n", rank, consistency_check, consistency_check_PD, node_count);
  // }
  PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "outputTXT"
PetscErrorCode outputTXT (PARAMETERS *par, ParticleManager &manager)
{


  PetscInt count=par->numNodes;
  PetscInt outCount=par->stepNumber;
  PetscInt counter=par->stepNumber/par->FreqResults;

  PetscInt num_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ostringstream convert;
  ostringstream convert1;
  ostringstream convert2;
  ostringstream convert3;
  ostringstream convert4;
  ostringstream convert5;
  PetscInt i,j;
  count = 0;
  int count1 = 0;
  BGL_FORALL_VERTICES(v,manager.graph,Graph){
    FieldData fd = get(fd_property,v);
    ParticleInfo info = get(info_property,v);
    if(!info.isTask){
        count++;
        // Geo << count1 <<  "  " << scientific << setprecision(4) << info.currentCoord[0] << "  " << info.currentCoord[1] << "  " << info.currentCoord[2] <$
        }
    }

  if (count > 0){
  if (PETSC_TRUE){

    // ##################################################
    //                  Geometry File
    // ##################################################

    convert << "Meshless." << rank << "." << counter << ".geo";
    string fname = convert.str();
    ofstream Geo;
    Geo.open(fname.c_str());

    Geo << "Meshless" << endl;
    Geo << "node" << endl;
    Geo << "node id given" << endl;
    Geo << "element id given" << endl;
    Geo << "coordinates" << endl;
    Geo << count << endl;

           BGL_FORALL_VERTICES(v,manager.graph,Graph){
           FieldData fd = get(fd_property,v);
           ParticleInfo info = get(info_property,v);
            if(!info.isTask){
                  count1++;
                  Geo << count1 <<  "  " << scientific << setprecision(4) << info.currentCoord[0] << "  " << info.currentCoord[1] << "  " << info.currentCoord[2] << endl;
            }
           }
    Geo << "part    1" << endl;
    Geo << "todo" << endl;
    Geo << "point" << endl;
    Geo << count << endl;
      for(j=1;j<=count;j++){
         Geo << j << " " << j << endl;
      }
      Geo.close();


      // ##################################################
      //                  Case File
      // ##################################################

      if (outCount == 0){
        convert1 << "Meshfree." << rank << ".case";
        fname = convert1.str();
        ofstream Case;
        Case.open(fname.c_str());

      PetscInt numSteps = (int) (par->finalTime/par->timeStep);

      Case << "#BOF: meshless.case" << endl;
      Case << endl;
      Case << "FORMAT" << endl;
      Case << endl;
      Case << "type: ensight" << endl;
      Case << endl;
      Case << "GEOMETRY" << endl;
      Case << endl;
      Case << "model: 1 Meshless." << rank << ".*.geo" << endl;
      Case << endl;
      Case << "VARIABLE" << endl;
      Case << endl;
      Case << "scalar per node: 1 Density Density." << rank << ".*.res" << endl;
      Case << "vector per node: 1 Velocity Velocity." << rank << ".*.res" << endl;
      Case << "scalar per node: 1 Pressure Pressure." << rank << ".*.res" << endl;
      Case << endl;
      Case << "TIME" << endl;
      Case << endl;
      Case << "time set: 1" << endl;
      Case << "number of steps: " << numSteps/par->FreqResults << endl;
      Case << "filename start number: 0" << endl;
      Case << "filename increment: " << "1" << endl;
      Case << "time values:" << endl;

      PetscInt counter1=0;
      for(i=0; i< (par->finalTime/par->timeStep);i++){
        if (i % par->FreqResults ==0){
          Case << counter1 << " ";
          counter1++;
        }
        if ((i+1) % (10*par->FreqResults) == 0) Case << endl;
      }
      Case.close();
      }


//  // ##################################################
//  //                  Density File
//  // ##################################################

      convert2 << "Density." << rank << "." << counter << ".res";
      fname = convert2.str();
      ofstream Density;
      Density.open(fname.c_str());

      Density << "Density" << endl;
      PetscInt counter2 = 0;
           BGL_FORALL_VERTICES(v,manager.graph,Graph){
             ParticleInfo info = get(info_property,v);
             FieldData fd = get(fd_property,v);
            if(!info.isTask){
              Density << scientific << setprecision(4) << fd.nodalDensity << "  ";
              counter2++;
              if (((counter2) % 6 == 0)) Density << endl;
            }
           }
           Density.close();

convert4 << "Pressure." << rank << "." << counter << ".res";
fname = convert4.str();
ofstream Pressure;
Pressure.open(fname.c_str());
Pressure << "Pressure" << endl;
counter2 = 0;
BGL_FORALL_VERTICES(v,manager.graph,Graph){
    ParticleInfo info = get(info_property,v);
    FieldData fd = get(fd_property,v);
    if(!info.isTask){
      Pressure << scientific << setprecision(4) << fd.nodalPressure << "  ";
      counter2++;
      if (((counter2) % 6 == 0)) Pressure << endl;
        }
      }
    Pressure.close();

//  // ##################################################
//  //                  Velocity File
//  // ##################################################

          convert5 << "Velocity." << rank << "." << counter << ".res";
          fname = convert5.str();
          ofstream Velocity;
          Velocity.open(fname.c_str());
          PetscInt counter3 = 0;
          Velocity << "Velocity" << endl;
               BGL_FORALL_VERTICES(v,manager.graph,Graph){
                 ParticleInfo info = get(info_property,v);
                 FieldData fd = get(fd_property,v);
                if(!info.isTask){
                  Velocity << scientific << setprecision(4) << fd.totalPhysicalVelocity[0] << "  " << fd.totalPhysicalVelocity[1] << "  " << fd.totalPhysicalVelocity[2] << "  "  ;
                  counter3++;
                  if (((counter3) % 2 == 0)) Velocity << endl;
                }

               }
               Velocity.close();

  }
}
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "OutputRestarts"
PetscErrorCode OutputRestarts(PARAMETERS *par,Vec U,Vec V,ParticleManager &manager)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i,j;
  PetscInt count=par->numNodes;


  PetscInt num_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ostringstream convert;
  ostringstream convert1;
  ostringstream convert2;
  ostringstream convert3;
  ostringstream convert4;
  ostringstream convert5;
  ostringstream convert6;
  ostringstream convert7;
  ostringstream convert8;
  ostringstream convert9;
  ostringstream convert10;
  ostringstream convert11;

  string fname;

  PetscInt counter = 0; //Number of nodes stored in current rank
  BGL_FORALL_VERTICES(v,manager.graph,Graph){
 	 ParticleInfo info = get(info_property,v);
 	 if(!info.isTask){
 		 counter++;
 	 }
  }

  // ##################################################
  //                  Velocity File
  // ##################################################

    MPI_Comm comm;
    PetscViewer viewer;
    ierr = PetscObjectGetComm((PetscObject)U,&comm);CHKERRQ(ierr);
    char filenameRestart[256];
    sprintf(filenameRestart,"RestartU%d.dat",par->stepNumber);
    ierr = PetscViewerBinaryOpen(comm,filenameRestart,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(U,viewer);CHKERRQ(ierr);

    ierr = PetscObjectGetComm((PetscObject)V,&comm);CHKERRQ(ierr);

    sprintf(filenameRestart,"RestartV%d.dat",par->stepNumber);
    ierr = PetscViewerBinaryOpen(comm,filenameRestart,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(V,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);

		convert5 << "RestartVelocity." << rank << "." << par->stepNumber << ".dat";
		fname = convert5.str();
		ofstream Velocity;
		Velocity.open(fname, ios::out | ios::binary);

		// ##################################################
		//                  Acceleration File
		// ##################################################

		  convert4 << "RestartAcceleration." << rank << "." << par->stepNumber << ".dat";
		  fname = convert4.str();
		  ofstream Acceleration;
		  Acceleration.open(fname, ios::out | ios::binary);

		// ##################################################
		//                  Stress File
		// ##################################################

		convert3 << "RestartStress." << rank << "." << par->stepNumber << ".dat";
		fname = convert3.str();
		ofstream Stress;
		Stress.open(fname, ios::out | ios::binary);

		// ##################################################
		//                  Strain File
		// ##################################################

		convert2 << "RestartStrain." << rank << "." << par->stepNumber << ".dat";
		fname = convert2.str();
		ofstream Strain;
		Strain.open(fname, ios::out | ios::binary);

		// ##################################################
		//                  Geometry File
		// ##################################################

		convert1 << "RestartGeo." << rank << "." << par->stepNumber << ".dat";
		fname = convert1.str();
		ofstream Geo;
		Geo.open(fname, ios::out | ios::binary);

		// ##################################################
		//                  Damage File
		// ##################################################

		convert << "RestartDamage." << rank << "." << par->stepNumber << ".dat";
		fname = convert.str();
		ofstream Damage;
		Damage.open(fname, ios::out | ios::binary);

		// ##################################################
		//               Damage Threshold File
		// ##################################################

		convert6 << "RestartThreshold." << rank << "." << par->stepNumber << ".dat";
		fname = convert6.str();
		ofstream Threshold;
		Threshold.open(fname, ios::out | ios::binary);

		// ##################################################
		//             Deformation Gradient File
		// ##################################################

		convert7 << "RestartDefGrad." << rank << "." << par->stepNumber << ".dat";
		fname = convert7.str();
		ofstream DefGrad;
		DefGrad.open(fname, ios::out | ios::binary);

		// ##################################################
		//                  Material File
		// ##################################################

		convert8 << "RestartMat." << rank << "." << par->stepNumber << ".dat";
		fname = convert8.str();
		ofstream Mat;
		Mat.open(fname, ios::out | ios::binary);

		// ##################################################
		//                  Volume File
		// ##################################################

		convert9 << "RestartVol." << rank << "." << par->stepNumber << ".dat";
		fname = convert9.str();
		ofstream Vol;
		Vol.open(fname, ios::out | ios::binary);

		// ##################################################
		//                  Boundary File
		// ##################################################

		convert10 << "RestartBound." << rank << "." << par->stepNumber << ".dat";
		fname = convert10.str();
		ofstream Bound;
		Bound.open(fname, ios::out | ios::binary);

    // ##################################################
    //              Number of particles File
    // ##################################################

    convert11 << "RestartNum." << rank << "." << par->stepNumber << ".dat";
    fname = convert11.str();
    ofstream Num;
    Num.open(fname, ios::out | ios::binary);
		Num.write( (char*)&counter, sizeof(int));


		if (counter > 0){
		   BGL_FORALL_VERTICES(v,manager.graph,Graph){
			 ParticleInfo info = get(info_property,v);
			 FieldData fd = get(fd_property,v);
			if(!info.isTask){

				Velocity.write( (char*)&fd.totalPhysicalVelocity[0], sizeof(double));
				Velocity.write( (char*)&fd.totalPhysicalVelocity[1], sizeof(double));
				Velocity.write( (char*)&fd.totalPhysicalVelocity[2], sizeof(double));

				Acceleration.write( (char*)&fd.totalPhysicalAcceleration[0], sizeof(double));
				Acceleration.write( (char*)&fd.totalPhysicalAcceleration[1], sizeof(double));
				Acceleration.write( (char*)&fd.totalPhysicalAcceleration[2], sizeof(double));

				Stress.write( (char*)&fd.totalStress0[0], sizeof(double));
				Stress.write( (char*)&fd.totalStress0[1], sizeof(double));
				Stress.write( (char*)&fd.totalStress0[2], sizeof(double));
				Stress.write( (char*)&fd.totalStress0[3], sizeof(double));
				Stress.write( (char*)&fd.totalStress0[4], sizeof(double));
				Stress.write( (char*)&fd.totalStress0[5], sizeof(double));

				Strain.write( (char*)&fd.totalStrain0[0], sizeof(double));
				Strain.write( (char*)&fd.totalStrain0[1], sizeof(double));
				Strain.write( (char*)&fd.totalStrain0[2], sizeof(double));
				Strain.write( (char*)&fd.totalStrain0[3], sizeof(double));
				Strain.write( (char*)&fd.totalStrain0[4], sizeof(double));
				Strain.write( (char*)&fd.totalStrain0[5], sizeof(double));

				Geo.write( (char*)&info.currentCoord[0], sizeof(double));
				Geo.write( (char*)&info.currentCoord[1], sizeof(double));
				Geo.write( (char*)&info.currentCoord[2], sizeof(double));

				DefGrad.write( (char*)&fd.currentDeformationGradient[0], sizeof(double));
				DefGrad.write( (char*)&fd.currentDeformationGradient[1], sizeof(double));
				DefGrad.write( (char*)&fd.currentDeformationGradient[2], sizeof(double));
				DefGrad.write( (char*)&fd.currentDeformationGradient[3], sizeof(double));
				DefGrad.write( (char*)&fd.currentDeformationGradient[4], sizeof(double));
				DefGrad.write( (char*)&fd.currentDeformationGradient[5], sizeof(double));
				DefGrad.write( (char*)&fd.currentDeformationGradient[6], sizeof(double));
				DefGrad.write( (char*)&fd.currentDeformationGradient[7], sizeof(double));
				DefGrad.write( (char*)&fd.currentDeformationGradient[8], sizeof(double));

				Damage.write( (char*)&fd.damage0, sizeof(double));

				Threshold.write( (char*)&fd.ductile_threshold0, sizeof(double));
				Threshold.write( (char*)&fd.brittle_threshold0, sizeof(double));

				Mat.write( (char*)&fd.material, sizeof(int));

        //Volume file includes information for volume update
        //[volume, volume_initial, nodalDensity, nodalDensityInitial]
				Vol.write( (char*)&fd.nodalVolume, sizeof(double));
        //PetscPrintf(PETSC_COMM_WORLD," Restart Nodal Vol %e \n", fd.nodalVolume);
        Vol.write( (char*)&fd.nodalVolumeInitial, sizeof(double));
        Vol.write( (char*)&fd.nodalDensity, sizeof(double));
        Vol.write( (char*)&fd.nodalDensityInitial, sizeof(double));

				Bound.write( (char*)&fd.Boundary, sizeof(int));
			}
		   }
    }
		   Velocity.close();
		   Acceleration.close();
		   Stress.close();
		   Strain.close();
		   Geo.close();
		   Damage.close();
		   Threshold.close();
		   DefGrad.close();
		   Mat.close();
		   Vol.close();
		   Bound.close();
		   Num.close();

  PetscFunctionReturn(0);
}

//// End I/O ////




//########## Main Function #############//

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[]) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx user;
  PetscInt i, j;
  PetscInt num_proc;
  PetscMPIInt rank;
  int mpi_id = 0;
  int mpi_size = 1;
  char version[128];

  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
  Teuchos::RCP<Epetra_Comm> epetraComm;

  #ifdef HAVE_MPI
    MPI_Comm_size(PETSC_COMM_WORLD, &num_proc);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_id);
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    epetraComm = Teuchos::RCP<Epetra_Comm>(new Epetra_MpiComm(PETSC_COMM_WORLD));
  #else
    epetraComm = Teuchos::RCP<Epetra_Comm>(new Epetra_SerialComm);
  #endif

  // Banner //
  PetscGetVersion(version,sizeof(version));
  if(mpi_id == 0){
    cout << "\n-- Meshfree Geometry Visualization --" << endl ;
    cout << "Petsc Version\n" << version<< endl ;
    cout << "Using Boost "
          << BOOST_VERSION / 100000     << "."
          << BOOST_VERSION / 100 % 1000 << "."
          << BOOST_VERSION % 100
          << endl;
    if(mpi_size > 1){
      cout << "MPI initialized on " << mpi_size << " processors.\n" << endl;}
  }

  PARAMETERS *par;
  ierr = PetscMalloc1(sizeof(*par),&par);CHKERRQ(ierr);

  //Material properties of Water at shallow depth
  user.mu       = 8.9e-4;
  user.lamda    = -2.0*user.mu/3.0;
  user.kappa    = 0.6;
  user.Cv       = 4180.0;
  user.Cp       = 4181.1;
  user.p0       = 100000.0;
  user.max_its  = 3;
  user.thickness = 0.00126;

  /* Set discretization options */
  PetscInt p=2, C=PETSC_DECIDE;

  // MPI initialization
  boost::mpi::environment env(argc,argv);
  // create an empty particle manager
  ParticleManager manager = ParticleManager();
  // seed random number generator based on MPI rank
  srand(manager.myRank);

  par->initialTime  	   = 0.0;
  par->finalTime    	   = 1.0e-3;
  par->timeStep	         = 2.5e-7;
  par->currentTime  	   = par->initialTime;
  par->gamma			       = 0.5;
  par->stepNumber        = 0;
  par->FreqResults       = 10;
  par->densityRDX        = 1770.0;
  par->density           = 1.42e3;//7.8e3; //Properties of composite plate
  par->youngModulus      = 78.4e9;//200e9;
  par->poissonRatio      = 0.039;//0.3;

  // Parameters for Penalty based approach
  par->penaltyConstant = 1.0;
  par->DamageModeling = PETSC_FALSE;
  par->damageCriticalStress = 8.0e10;
  par->damageCriticalEpsilonPlastic = 0.2;
  par->thresholdDamageForPenalty = 0.9;

  user.TimeRestart       = 0;
  user.StepRestart       = 0;
  user.FreqRestarts      = 10000;

  user.spacing  = 0.0000000000000000;
  user.Lx       = 1.8+user.spacing;
  user.Ly       = 1.8+user.spacing;
  user.Lz       = 1.8+user.spacing;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","PhaseFieldCrystal2D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-StepRestart","Step of the initial solution from file",__FILE__,user.StepRestart,&user.StepRestart,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (C == PETSC_DECIDE) C = p-1;

  //IGA
  PetscInt dim = 3;
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,dim);CHKERRQ(ierr);
  ierr = IGASetDof(iga,5);CHKERRQ(ierr);
  ierr = IGAGetDim(iga,&dim);CHKERRQ(ierr);
  ierr = IGARead(iga,"./Geometry.dat");CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  ierr = IGAWrite(iga,"igaF.dat");CHKERRQ(ierr);
  user.iga = iga;

  //Create solution vector V (velocities) and A (accelerations) and dispacements D
  PetscReal t=0;
  ierr = IGACreateVec(iga,&user.V0);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.A0);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.D0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.V0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.A0);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.D0);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&user.dA);CHKERRQ(ierr);
  ierr = VecZeroEntries(user.dA);CHKERRQ(ierr);

  ierr = PetscMalloc1(iga->proc_sizes[0],&user.processor_numElX);CHKERRQ(ierr);
  ierr = PetscMalloc1(iga->proc_sizes[1],&user.processor_numElY);CHKERRQ(ierr);
  ierr = PetscMalloc1(iga->proc_sizes[2],&user.processor_numElZ);CHKERRQ(ierr);

  PetscInt tempX[iga->proc_sizes[0]];
  PetscInt tempY[iga->proc_sizes[1]];
  PetscInt tempZ[iga->proc_sizes[2]];

    for (j=0;j<iga->proc_sizes[0];j++){
      user.processor_numElX[j] = 0;
      tempX[j] = 0;
    }
    for (j=0;j<iga->proc_sizes[1];j++){
      user.processor_numElY[j] = 0;
      tempY[j] = 0;
    }
    for (j=0;j<iga->proc_sizes[2];j++){
      user.processor_numElZ[j] = 0;
      tempZ[j] = 0;
    }

    tempX[iga->proc_ranks[0]] = iga->elem_width[0];
    tempY[iga->proc_ranks[1]] = iga->elem_width[1];
    tempZ[iga->proc_ranks[2]] = iga->elem_width[2];

  ierr = MPI_Allreduce(&tempX,user.processor_numElX,iga->proc_sizes[0],MPI_INT,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&tempY,user.processor_numElY,iga->proc_sizes[1],MPI_INT,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&tempZ,user.processor_numElZ,iga->proc_sizes[2],MPI_INT,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);

  for (j=0;j<iga->proc_sizes[0];j++) user.processor_numElX[j] /= iga->proc_sizes[1]*iga->proc_sizes[2];
  for (j=0;j<iga->proc_sizes[1];j++) user.processor_numElY[j] /= iga->proc_sizes[0]*iga->proc_sizes[2];
  for (j=0;j<iga->proc_sizes[2];j++) user.processor_numElZ[j] /= iga->proc_sizes[0]*iga->proc_sizes[1];


  user.H[0] = user.Lx/iga->axis[0]->nel;
  user.H[1] = user.Ly/iga->axis[1]->nel;
  user.H[2] = user.Lz/iga->axis[2]->nel;


    par->stepNumber  = 0;
    ierr = input(par,manager,&user);CHKERRQ(ierr);
    manager.sync();
    MPI_Allreduce(&user.totalInitialExplosiveVolume, &totalInitialExplosiveVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&par->numNodes, &totalNumNodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Read %d Immersed Particles from %d Foreground files\n", totalNumNodes, mpi_size);
    manager.connectVertsToTasks(false, &user);
    manager.sync();


    ierr = FormInitialCondition(iga,t,user.V0,&user);CHKERRQ(ierr);
     // Dump Initial Solution
     char filename[256];
     sprintf(filename,"velS%d.dat",par->stepNumber);
     ierr = IGAWriteVec(user.iga,user.V0,filename);CHKERRQ(ierr);

  if (par->stepNumber == 0){
    BGL_FORALL_VERTICES(v,manager.graph,Graph){
      FieldData fd = get(fd_property,v);
      ParticleInfo info = get(info_property,v);
      if(!info.isTask){

        fd.Inside = 1;
        fd.Boundary = 0;
        fd.damage = 0;
        fd.damage0 = 0;

      for (j=0;j<3;j++){
        fd.totalPhysicalAcceleration[j] = 0.0;
        fd.totalPhysicalDisplacement[j] = 0.0;
        fd.totalPhysicalVelocity[j] = 0.0;

        fd.totalPhysicalDisplacementOldStep[j] = 0.0;
        fd.totalPhysicalVelocityOldStep[j] = 0.0;
        fd.totalPhysicalAccelerationOldStep[j] = 0.0;

        fd.totalPhysicalDisplacementOldIteration[j] = 0.0;
        fd.totalPhysicalVelocityOldIteration[j] = 0.0;
        fd.totalPhysicalAccelerationOldIteration[j] = 0.0;

        fd.AccelerationIncrement[j] = 0.0;
        fd.interpolatedVelocity[j] = 0.0;

        fd.referenceCoord[j] = fd.currentCoord[j];

      }

      for(j = 0 ; j < 9 ; j++){
      fd.currentDeformationGradient[j] = 0.0;
      fd.alphaDeformationGradient[j] = 0.0;
      fd.DeformationGradientOld[j] = 0.0;

      fd.velocityGradient[j] = 0.0;
      }

      fd.DeformationGradientOld[0] = 1.0;
      fd.DeformationGradientOld[4] = 1.0;
      fd.DeformationGradientOld[8] = 1.0;

      fd.alphaDeformationGradient[0] = 1.0;
      fd.alphaDeformationGradient[4] = 1.0;
      fd.alphaDeformationGradient[8] = 1.0;

      fd.currentDeformationGradient[0] = 1.0;
      fd.currentDeformationGradient[4] = 1.0;
      fd.currentDeformationGradient[8] = 1.0;

      fd.determinantAlphaDeformationGradient = 1.0;
      fd.determinantCurrentDeformationGradient = 1.0;
      fd.nodalDensity = 1.0;
      fd.nodalPressure = 1.0;

      put(fd_property,v,fd);
    }
  }
}
manager.sync();

  const PetscScalar *arrayU;
  const PetscScalar *arrayV;
  Vec               localV;
  Vec               localU;

    PetscPrintf(PETSC_COMM_WORLD,"######################################################### \n");
    PetscPrintf(PETSC_COMM_WORLD," Writing Meshfree Geometry files... \n");
    PetscPrintf(PETSC_COMM_WORLD,"######################################################### \n");

    par->currentTime+=par->timeStep;
    user.stepNumber  =par->stepNumber;

    ierr = outputTXT(par,manager);CHKERRQ(ierr);

PetscPrintf(PETSC_COMM_WORLD, "Complete\n");
#ifdef HAVE_MPI
  PetscFree(par);CHKERRQ(ierr);
  IGADestroy(&iga);CHKERRQ(ierr);
  PetscFinalize();CHKERRQ(ierr);
#endif
PetscFunctionReturn(0);
}
