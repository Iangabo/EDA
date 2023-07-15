#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <limits>
#include <chrono>

// En esta clase vamos a recibir las claves y distancias para luego pasarlas a arbol RDB
class DataObject {
public:
    unsigned long long int hilbert_key;
    std::vector<double> distances; // Vector de distancias
    std::vector<double> coordinates; // Vector de coordenadas

    DataObject() : hilbert_key(0) {} // Constructor predeterminado

    DataObject(unsigned long long int key, double dist) : hilbert_key(key) {
        // Inicializa distances con dist
        distances.push_back(dist);
    }

    DataObject(unsigned long long int key, const std::vector<double>& dist)
        : hilbert_key(key), distances(dist) {}

    DataObject(unsigned long long int key, const std::vector<double>& coords, double dist)
        : hilbert_key(key), coordinates(coords), distances(dist) {}

    bool operator==(const DataObject& other) const {
        // Compara las coordenadas de los objetos
        return coordinates == other.coordinates;
    }
};


//Para puntos n dimensionales
struct PointND {
    std::vector<double> dimensions;
    std::string songName;
};

//Para puntos 2D que generara partitionPoints
struct Point2D {
    double x;
    double y;
};

//Funcion para leer los puntos desde el archivo CSV
std::vector<PointND> readPointsFromCSV(const std::string& filename) {
    std::vector<PointND> points;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cout << "Error al abrir el archivo." << std::endl;
        return points;
    }

    std::string line;
    std::getline(file, line); // Ignorar la primera línea

    int lineNumber = 1;
    while (std::getline(file, line)) {
        PointND point;
        std::istringstream iss(line);
        std::string value;

        std::vector<double> coordinates; // Nuevo vector para almacenar las coordenadas

        while (std::getline(iss, value, '\t')) {
            try {
                double dimension = std::stod(value);
                point.dimensions.push_back(dimension);
                coordinates.push_back(dimension); // Agrega la coordenada al vector
            }
            catch (const std::exception& e) {
                // Si no se puede convertir a double, se considera el último valor como el songName
                point.songName = value;
            }
        }

        // Crea un objeto DataObject con las coordenadas y las distancias
        DataObject data(0, coordinates, 0.0); // Actualiza los valores apropiados

        points.push_back(point);
        //std::cout << "Punto leido: " << lineNumber << std::endl;
        lineNumber++;
    }

    file.close();
    return points;
}
//Partimos un punto n dimensional en punto sde 2D para sacar sus claves de hilbert
std::vector<std::vector<Point2D>> partitionPoints(const std::vector<PointND>& points, int dimensionsPerPartition) {
    std::vector<std::vector<Point2D>> partitions;
    int totalDimensions = points[0].dimensions.size();

    // Crear particiones
    for (int d = 0; d < totalDimensions; d += dimensionsPerPartition) {
        std::vector<Point2D> partition;
        for (const PointND& point : points) {
            Point2D point2D;
            point2D.x = point.dimensions[d];
            point2D.y = (d + 1 < totalDimensions) ? point.dimensions[d + 1] : 0;
            partition.push_back(point2D);
        }
        partitions.push_back(partition);
    }

    return partitions;
}

// Calculamos la distancia euclidiana
double calculateEuclideanDistance(const PointND& queryPoint, const DataObject& object) {
    const std::vector<double>& queryCoordinates = queryPoint.dimensions;
    const std::vector<double>& objectCoordinates = object.coordinates;

    if (queryCoordinates.size() != objectCoordinates.size()) {
        // Manejar un error o una situación inválida
        return 0.0;
    }

    double distance = 0.0;
    for (size_t i = 0; i < queryCoordinates.size(); ++i) {
        double diff = queryCoordinates[i] - objectCoordinates[i];
        distance += diff * diff;
    }

    return std::sqrt(distance);
}


//Con esta seccion vamos a calcular las claves de hilbert

struct cuadranteHilbert {
    int primer = 0;
    int segundo = 1;
    int tercer = 2;
    int cuarto = 3;
};

unsigned long long int hilbertKey(Point2D p, int n, double xMin, double xMax, double yMin, double yMax, cuadranteHilbert cuadrante) {
    if (n <= 0)
        return 0;

    long double xMid = (xMin + xMax) / 2.0;
    long double yMid = (yMin + yMax) / 2.0;

    unsigned int cuad;
    unsigned int real;
    if (p.x < xMid)
    {
        if (p.y < yMid)
        {
            cuad = cuadrante.primer;
            real = 0;
        }
        else
        {
            cuad = cuadrante.segundo;
            real = 1;
        }

    }
    else
    {
        if (p.y < yMid)
        {
            cuad = cuadrante.cuarto;
            real = 3;
        }
        else
        {
            cuad = cuadrante.tercer;
            real = 2;
        }

    }

    cuadranteHilbert aux;

    if (cuad == 0) {
        int tmp = aux.segundo;
        aux.segundo = aux.cuarto;
        aux.cuarto = tmp;
    }

    if (cuad == 3) {
        int tmp = aux.primer;
        aux.primer = aux.tercer;
        aux.tercer = tmp;
    }

    double subSquareXMin, subSquareXMax, subSquareYMin, subSquareYMax;
    switch (real) {
    case 0:
        subSquareXMin = xMin;
        subSquareXMax = xMid;
        subSquareYMin = yMin;
        subSquareYMax = yMid;
        break;
    case 1:
        subSquareXMin = xMin;
        subSquareXMax = xMid;
        subSquareYMin = yMid;
        subSquareYMax = yMax;
        break;
    case 2:
        subSquareXMin = xMid;
        subSquareXMax = xMax;
        subSquareYMin = yMid;
        subSquareYMax = yMax;
        break;
    default:
        subSquareXMin = xMid;
        subSquareXMax = xMax;
        subSquareYMin = yMin;
        subSquareYMax = yMid;
        break;
    }

    unsigned long long int subSquareKey = hilbertKey(p, n - 1, subSquareXMin, subSquareXMax, subSquareYMin, subSquareYMax, aux);
    unsigned long long int key;
    if (n == 1)
    {
        if (cuad == 0)
        {
            key = 0 + subSquareKey;
        }
        else if (cuad == 1) {
            key = 1 + subSquareKey;
        }
        else if (cuad == 2) {
            key = 2 + subSquareKey;
        }
        else if (cuad == 3) {
            key = 3 + subSquareKey;
        }
    }
    else
    {
        if (cuad == 0) {
            key = 0 * (unsigned long long int)std::pow(2, (unsigned long long int)pow(2, n - 1)) + subSquareKey;
        }
        else if (cuad == 1) {
            key = 1 * (unsigned long long int)std::pow(2, (unsigned long long int)pow(2, n - 1)) + subSquareKey;
        }
        else if (cuad == 2) {
            key = 2 * (unsigned long long int)std::pow(2, (unsigned long long int)pow(2, n - 1)) + subSquareKey;
        }
        else if (cuad == 3) {
            key = 3 * (unsigned long long int)std::pow(2, (unsigned long long int)pow(2, n - 1)) + subSquareKey;
        }
    }
    return key;
}

//Con esta seccion vamos a construir los arboles RDB

//Nodo del arbol RDB
class NodoRDBTree {
public:
    bool isLeaf;
    std::vector<DataObject> objects;
    std::vector<NodoRDBTree*> children;
    std::vector<double> distances;

    NodoRDBTree(bool leaf = false) : isLeaf(leaf) {}
};
//Arbol RDB
class RDBTree {
private:
    int order; //orden del arbol RDB
    std::vector<DataObject> referenceobjectos;

public:
    NodoRDBTree* root;
    RDBTree() {
        order = 20;//Definimos el orden
        root = nullptr;
    }
    // Insertamos un objeto de datos en el arbol RDB
    void insert(const DataObject& data, const std::vector<double>& distances) {
        if (root == nullptr) {
            root = new NodoRDBTree(true);
            root->objects.push_back(data);
            if (root->isLeaf) {
                root->distances.push_back(distances.back());  // Almacenar la distancia correspondiente
            }
        }
        else {
            if (root->objects.size() == (order + 1)) {
                NodoRDBTree* newRoot = new NodoRDBTree(false);
                newRoot->children.push_back(root);
                splitChild(newRoot, 0, root, data);
                insertarNodoNoLleno(newRoot, data, distances);
                root = newRoot;
            }
            else {
                insertarNodoNoLleno(root, data, distances);
            }
        }
    }

    // Busqueda de una clave en el arbol
    NodoRDBTree* buscarKey(int hilbertKey) {
        if (root == nullptr) {
            std::cout << "Arbol vacio.\n";
            return nullptr;
        }
        else {
            return buscarKeyAux(root, hilbertKey);
        }
    }

    //Para verificar si es que se estan insertando las claves de hilbert en lo arboles RDB
    void mostrar() {
        if (root == nullptr) {
            std::cout << "Arbol vacio.\n";
        }
        else {
            mostrarAux(root);
        }
    }

    // Obtener una lista de objetos de datos del árbol RDB
    std::vector<DataObject> getObjects() {
        std::vector<DataObject> objects;
        getObjectsAux(root, objects);
        return objects;
    }

private:

    void splitChild(NodoRDBTree* parent, int index, NodoRDBTree* child, const DataObject& data) {
        // Creamos un nuevo nodo
        NodoRDBTree* newNode = new NodoRDBTree(child->isLeaf);
        if (child->isLeaf)
            // Si el nodo hijo es una hoja, insertamos el último objeto en el padre
            parent->objects.insert(parent->objects.begin() + index, child->objects[order - 1]);
        else
            // Si el nodo hijo no es una hoja, insertamos el objeto en el padre
            parent->objects.insert(parent->objects.begin() + index, child->objects[order]);

        // Si el nodo hijo no es una hoja, transferimos los nodos hijos restantes al nuevo nodo
        if (!child->isLeaf) {
            for (int i = order + 1; i < child->children.size(); ++i) {
                // Movemos los nodos hijos al nuevo nodo
                newNode->children.push_back(child->children[i]);
            }
            // Eliminamos los nodos hijos transferidos del nodo hijo original
            child->children.erase(child->children.begin() + order + 1, child->children.end());
        }

        //Insertamos el nuevo nodo como hijo del padre
        parent->children.push_back(newNode);

        if (!child->isLeaf)
        {
            // Si el nodo hijo no es una hoja, transferimos los objetos restantes al nuevo nodo
            for (int i = order + 1; i < child->objects.size(); ++i) {
                // Movemos los objetos al nuevo nodo
                newNode->objects.push_back(child->objects[i]);
            }
            // Eliminamos los objetos transferidos del nodo hijo original
            child->objects.erase(child->objects.begin() + order, child->objects.end());
        }
    }
    // Insertar en un nodo no lleno
    void insertarNodoNoLleno(NodoRDBTree* node, const DataObject& data, const std::vector<double>& distances) {
        int i = node->objects.size() - 1;
        if (node->isLeaf) {
            node->objects.push_back(DataObject(0, 0.0));
            node->objects[i + 1] = data;
            node->distances.push_back(distances.back());  // Almacenar la distancia correspondiente
        }
        else {
            i++;
            if (node->children[i]->objects.size() == order + 1) {
                splitChild(node, i, node->children[i], data);
                if (data.hilbert_key > node->objects[i].hilbert_key) {
                    i++;
                }
            }
            insertarNodoNoLleno(node->children[i], data, distances);

            if (node->children[i]->isLeaf == false && node->children[i]->children[0]->objects.size() == order + 1) {
                splitChild(node->children[i], 0, node->children[i]->children[0], data);
            }
        }
    }

    NodoRDBTree* buscarKeyAux(NodoRDBTree* node, int hilbertKey) {
        int i = 0;
        while (i < node->objects.size() && hilbertKey > node->objects[i].hilbert_key) {
            i++;
        }
        if (!node->children[i]->isLeaf)
            return buscarKeyAux(node->children[i], hilbertKey);
        else {
            return node->children[i];
        }
    }
    void mostrarAux(NodoRDBTree* node) {
        if (node->isLeaf) {
            std::cout << "Hoja: ";
            std::vector<DataObject> sortedObjects = node->objects;
            std::sort(sortedObjects.begin(), sortedObjects.end(), [](const DataObject& obj1, const DataObject& obj2) {
                return obj1.hilbert_key < obj2.hilbert_key;
                });
            for (size_t i = 0; i < sortedObjects.size(); ++i) {
                std::cout << "Clave: " << sortedObjects[i].hilbert_key;
                // Imprimir las distancias almacenadas en el nodo actual
                std::cout << ", Distancias: ";
                for (const auto& distance : sortedObjects[i].distances) {
                    std::cout << distance << " ";
                }

                std::cout << std::endl;
            }
        }
        else {
            std::cout << "Rama: ";
            std::vector<DataObject> sortedObjects = node->objects;
            std::sort(sortedObjects.begin(), sortedObjects.end(), [](const DataObject& obj1, const DataObject& obj2) {
                return obj1.hilbert_key < obj2.hilbert_key;
                });
            for (const auto& obj : sortedObjects) {
                std::cout << obj.hilbert_key << " ";
            }
            std::cout << std::endl;

            for (const auto& child : node->children) {
                mostrarAux(child);
            }
        }
    }

    // Función auxiliar para realizar el recorrido en orden y obtener los objetos
    void getObjectsAux(NodoRDBTree* node, std::vector<DataObject>& objects) {
        if (node == nullptr) {
            return;
        }

        if (node->isLeaf) {
            for (const auto& object : node->objects) {
                objects.push_back(object);
            }
        }
        else {
            for (const auto& child : node->children) {
                getObjectsAux(child, objects);
            }
        }
    }

};

// Clase HDIndex
class HDIndex {
public:
    RDBTree rdbTree;
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::lowest();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::lowest();
    std::string csvFilename;
    std::vector<PointND> points;
    std::vector<PointND> refPoints;
    std::vector<std::vector<double>> distances;
    std::vector<unsigned long long int> hilbertKeys;

    HDIndex(const std::string& filename) : csvFilename(filename) {}

    void buildIndex() {

        // Leer los puntos desde el archivo CSV
        std::vector<PointND> points = readPointsFromCSV(csvFilename);

        // Calcular los valores de xMin, xMax, yMin, yMax a partir de los puntos
        for (const auto& point : points) {
            for (const auto& dimension : point.dimensions) {
                if (dimension < xMin) {
                    xMin = dimension;
                }
                if (dimension > xMax) {
                    xMax = dimension;
                }
            }
            if (point.dimensions.size() > 1) {
                double y = point.dimensions[1];
                if (y < yMin) {
                    yMin = y;
                }
                if (y > yMax) {
                    yMax = y;
                }
            }
        }

        // Particionar los puntos
        std::vector<std::vector<Point2D>> partitions = partitionPoints(points, 2);

        // Calcular las claves de Hilbert
        std::vector<unsigned long long int> hilbertKeys;
        for (const auto& partition : partitions) {
            for (const auto& p : partition) {
                unsigned long long int key = hilbertKey(p, 6, xMin, xMax, yMin, yMax, cuadranteHilbert());
                hilbertKeys.push_back(key);
            }
        }

        // Crear objetos DataObject con las claves de Hilbert y el songName
        std::vector<DataObject> dataObjects;
        for (size_t i = 0; i < hilbertKeys.size(); ++i) {
            DataObject data(hilbertKeys[i], 0.0);
            dataObjects.push_back(data);
        }

        // Estimar la distancia máxima entre cualquier par de puntos en el conjunto de datos
        int r = rand() % points.size();
        int j = 0;
        double max_dist = 0;
        for (int t = 0; t < 10; t++) {
            max_dist = 0;
            for (int k = 0; k < points.size(); k++) {
                double dist = 0;
                for (int d = 0; d < points[k].dimensions.size(); d++) {
                    dist += (points[r].dimensions[d] - points[k].dimensions[d]) * (points[r].dimensions[d] - points[k].dimensions[d]);
                }
                if (dist > max_dist) {
                    max_dist = dist;
                    j = k;
                }
            }
            r = j;
        }
        max_dist = sqrt(max_dist);

        // Seleccionar puntos de referencia utilizando el algoritmo SSS
        refPoints.clear();
        refPoints.push_back(points[r]);
        double f = 0.8;
        int m = 3; // Definir el valor de m según tus necesidades

        while (refPoints.size() < m) {
            for (int i = 0; i < points.size(); i++) {
                bool is_valid = true;
                for (int j = 0; j < refPoints.size(); j++) {
                    double dist = 0;
                    for (int d = 0; d < points[i].dimensions.size(); d++) {
                        dist += (points[i].dimensions[d] - refPoints[j].dimensions[d]) * (points[i].dimensions[d] - refPoints[j].dimensions[d]);
                    }
                    dist = sqrt(dist);
                    if (dist < f * max_dist) {
                        is_valid = false;
                        break;
                    }
                }
                if (is_valid) {
                    refPoints.push_back(points[i]);
                    break;
                }
            }
            f *= 0.95;
        }

        // Calcular las distancias a los objetos de referencia para cada punto
        std::vector<std::vector<double>> distances;
        for (const auto& point : points) {
            std::vector<double> pointDistances;
            for (const auto& refPoint : refPoints) {
                double dist = 0;
                for (int d = 0; d < point.dimensions.size(); d++) {
                    dist += (point.dimensions[d] - refPoint.dimensions[d]) * (point.dimensions[d] - refPoint.dimensions[d]);
                }
                pointDistances.push_back(sqrt(dist));
            }
            distances.push_back(pointDistances);
        }

        auto start = std::chrono::steady_clock::now();
        // Insertar en el arbol RDB
        for (size_t i = 0; i < dataObjects.size(); ++i) {
            const auto& currentDistances = distances[i]; // Obtener el vector de distancias para el punto actual
            this->rdbTree.insert(DataObject(dataObjects[i].hilbert_key, currentDistances), currentDistances);
        }

        //Imprimir el contenido del árbol RDB
        //rdbTree.mostrar();
        //medir e imprimir tiempos
        auto end = std::chrono::steady_clock::now();
        int ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Tiempo de indexacion de objetos en el arbol:" << "\n" << std::endl;

        std::cout << "Nanosegundos: "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
            << " ns \n" << std::endl;

        std::cout << "Microsegundos: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << " us \n" << std::endl;

        std::cout << "Milisegundos: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << " ms \n" << std::endl;

        std::cout << "Segundos: "
            << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
            << " sec \n" << std::endl;
    }
};

int main() {
    HDIndex hdIndex("songs_final.csv");
    hdIndex.buildIndex();

    // Punto de consulta
    PointND queryPoint;
    queryPoint.dimensions.push_back(0.85);
    queryPoint.dimensions.push_back(893);
    queryPoint.dimensions.push_back(5);
    queryPoint.dimensions.push_back(-4.783);
    queryPoint.dimensions.push_back(1);
    queryPoint.dimensions.push_back(623);
    queryPoint.dimensions.push_back(138);
    queryPoint.dimensions.push_back(4.14e-06);
    queryPoint.dimensions.push_back(3720000000000000);
    queryPoint.dimensions.push_back(391);
    queryPoint.dimensions.push_back(218.05);
    queryPoint.dimensions.push_back(98821);
    queryPoint.dimensions.push_back(4);
    queryPoint.songName = "Symbiote";

    // Obtener la lista de objetos del árbol RDB
    std::vector<DataObject> objects = hdIndex.rdbTree.getObjects();
    int numeroObjetos = objects.size();
    //std::cout << "Numero de objetos en el arbol RDB: " << numeroObjetos << std::endl;

    // Crear un nuevo vector para almacenar los candidatos refinados
    std::vector<DataObject> refinedCandidates;

    // Iterar sobre los objetos del vector objects
    for (const auto& object : objects) {
        double maxTriangularBound = 0.0;
        double maxPtolemaicBound = 0.0;

        // Calcular los límites inferiores utilizando la desigualdad triangular
        for (size_t i = 0; i < hdIndex.refPoints.size(); ++i) {
            double distQRef = calculateEuclideanDistance(queryPoint, objects[i]);
            double distObjectRef = calculateEuclideanDistance(PointND{ object.coordinates, "" }, objects[i]);
            double triangularBound = std::abs(distQRef - distObjectRef);
            if (triangularBound > maxTriangularBound) {
                maxTriangularBound = triangularBound;
            }
        }


        // Calcular los límites inferiores utilizando la desigualdad ptolemaica
        for (size_t i = 0; i < hdIndex.refPoints.size(); ++i) {
            for (size_t j = i + 1; j < hdIndex.refPoints.size(); ++j) {
                double distQRefI = calculateEuclideanDistance(queryPoint, objects[i]);
                double distObjectRefJ = calculateEuclideanDistance(PointND{ object.coordinates, "" }, objects[j]);
                double distQRefJ = calculateEuclideanDistance(queryPoint, objects[j]);
                double distObjectRefI = calculateEuclideanDistance(PointND{ object.coordinates, "" }, objects[i]);
                double distRefIRefJ = calculateEuclideanDistance(hdIndex.refPoints[i], objects[j]);
                double ptolemaicBound = std::abs((distQRefI * distObjectRefJ) - (distQRefJ * distObjectRefI)) / distRefIRefJ;
                if (ptolemaicBound > maxPtolemaicBound) {
                    maxPtolemaicBound = ptolemaicBound;
                }
            }
        }

        // Comprobar si el objeto cumple con los límites inferiores
        if (maxTriangularBound < object.distances.back() && maxPtolemaicBound < object.distances.back()) {
            refinedCandidates.push_back(object);
        }
    }

    /*
    // Imprimir los candidatos refinados
    std::cout << "Numero de candidatos refinados: " << refinedCandidates.size() << std::endl;
    for (const auto& refinedCandidate : refinedCandidates) {
        std::cout << "Clave de Hilbert: " << refinedCandidate.hilbert_key << std::endl;
        std::cout << "Distancias: ";
        for (const auto& distance : refinedCandidate.distances) {
            std::cout << distance << " ";
        }
        std::cout << std::endl;
    }
    */

    // Definir el valor de k (número de vecinos más cercanos a obtener)
    int k;
    std::cout << "Ingrese numero de vecinos cercanos que desea buscar en la Estructura HD-Index: ";
    std::cin >> k;


    auto start = std::chrono::steady_clock::now();
    // Ordenar los candidatos refinados de menor a mayor distancia euclidiana
    std::partial_sort(refinedCandidates.begin(), refinedCandidates.begin() + k, refinedCandidates.end(), [&](const DataObject& obj1, const DataObject& obj2) {
        double dist1 = calculateEuclideanDistance(queryPoint, obj1);
        double dist2 = calculateEuclideanDistance(queryPoint, obj2);
        return dist1 < dist2;
        });
    auto end = std::chrono::steady_clock::now();

    // Imprimir los primeros k vecinos más cercanos
    std::cout << "Los " << k << " vecinos mas cercanos son: " << "\n" << std::endl;
    for (int i = 0; i < k; ++i) {
        const DataObject& neighbor = refinedCandidates[i];
        std::cout << "Clave de Hilbert: " << neighbor.hilbert_key << std::endl;
        std::cout << "Distancias: ";
        for (const auto& distance : neighbor.distances) {
            std::cout << distance << " ";
        }
        std::cout << "\n" << std::endl;
    }

    int ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tiempo de busqueda de los k vecinos en la estructura HD-Index:" << "\n" << std::endl;

    std::cout << "Nanosegundos: "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << " ns \n" << std::endl;

    std::cout << "Microsegundos: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " us \n" << std::endl;

    std::cout << "Milisegundos: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms \n" << std::endl;

    std::cout << "Segundos: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " sec \n" << std::endl;

    return 0;
}