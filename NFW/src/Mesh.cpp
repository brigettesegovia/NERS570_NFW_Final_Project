#include "Mesh.hpp"

#include <stddef.h>
#include <vector>
#include <array>
#include <cassert>
#include <iostream>

using namespace std;

namespace NBody
{
    using namespace std;

    // constructors

    /**
     * @brief Default constructor for Mesh
     *
     * Initializes a mesh object.
     *
     * Defaulted in Mesh.hpp
     */

    /**
     * @brief Copy constructor for Mesh
     *
     * Creates a copy of an existing mesh.
     *
     * @param other: existing Mesh object to copy
     *
     * Defaulted in Mesh.hpp
     */

    /**
     * @brief Copy assignment operator for Mesh
     *
     * Copies the contents of an existing mesh into this mesh.
     *
     * @param other: existing Mesh object to copy
     * @return Mesh&: reference to the updated mesh
     *
     * Defaulted in Mesh.hpp
     */

    /**
     * @brief Move constructor for Mesh
     *
     * Transfers the contents of an existing mesh into a new mesh.
     *
     * @param other: existing Mesh object to move
     *
     * Defaulted in Mesh.hpp
     */

    /**
     * @brief Move assignment operator for Mesh
     *
     * Transfers the contents of an existing mesh into this mesh.
     *
     * @param other: existing Mesh object to move
     * @return Mesh&: reference to the updated mesh
     *
     * Defaulted in Mesh.hpp
     */

    // destructor

    /**
     * @brief Virtual destructor for Mesh
     *
     * Allows derived mesh classes to be safely destroyed
     * through a base class pointer.
     *
     * Defaulted in Mesh.hpp
     */

    // getters

    /**
     * @brief Get the position array of the mesh
     *
     * @return size_t***: pointer to the 3D mesh position array
     *
     * Implemented directly in Mesh.hpp
     */

    // abstract methods

    /**
     * @brief Calculate the gravitational potential across the mesh
     *
     * Pure virtual method.
     * Implementation provided by MeshSerial or MeshParallel.
     */

    /**
     * @brief Calculate the density across the mesh
     *
     * Pure virtual method.
     * Implementation provided by MeshSerial or MeshParallel.
     */

    /**
     * @brief Update the positions of particles contained in the mesh
     *
     * Pure virtual method.
     * Implementation provided by MeshSerial or MeshParallel.
     */

    /**
     * @brief Update the velocities of particles contained in the mesh
     *
     * Pure virtual method.
     * Implementation provided by MeshSerial or MeshParallel.
     */

    /**
     * @brief Initialize particles within the mesh
     *
     * Pure virtual method.
     * Implementation provided by MeshSerial or MeshParallel.
     */

}
