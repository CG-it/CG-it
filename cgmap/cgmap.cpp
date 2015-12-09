#include <tcl.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

/**
 *
 * +----------+
 * | CGit     |
 * +----------+
 *
 * CGit a VMD package to simplify creating coarse grained SDK
 * topologies.
 *
 * Copyright (c) 2013 by Chris MacDermaid <chris.macdermaid@gmail.com>
 * and Giacomo Fiorin <giacomo.fiorin@gmail.com>
 *
 *(Shinoda) Shinoda, DeVane, Klein, Mol Sim, 33, 27 (2007).
 *(DeVane) Shinoda, DeVane, Klein, Soft Matter, 4, 2453-2462 (2008).
 *
 * Some of the work here is used verbatim from Jerome Henin's clever
 * implementation of qwrap: https://github.com/jhenin/qwrap
 **/

extern "C" {
    int Cgmap_Init(Tcl_Interp *interp);
}

static int obj_Cgmap(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[]);

int parse_vector (Tcl_Obj * const obj, std::vector<float> &vec, Tcl_Interp *interp)
{
    Tcl_Obj **data;
    int num;
    double d;

    if (Tcl_ListObjGetElements(interp, obj, &num, &data) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing arguments", TCL_STATIC);
        return -1;
    }

    vec.resize(num);

    for (int i = 0; i < num; i++) {
        if (Tcl_GetDoubleFromObj(interp, data[i], &d) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Cgmap: error parsing vector element as floating-point", TCL_STATIC);
            return -1;
        }
        // Tcl gives us doubles, make them float
        vec[i] = float (d);
    }
    return num;
}

int parse_ivector (Tcl_Obj * const obj, std::vector<int> &vec, Tcl_Interp *interp, bool fromDouble)
{
    Tcl_Obj **data;
    int num;
    double d;

    if (Tcl_ListObjGetElements(interp, obj, &num, &data) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing arguments", TCL_STATIC);
        return -1;
    }

    vec.resize(num);

    if (fromDouble == false) {
        for (int i = 0; i < num; i++) {
            if (Tcl_GetIntFromObj(interp, data[i], &vec[i]) != TCL_OK) {
                Tcl_SetResult(interp, (char *) "Cgmap: error parsing vector element as integer", TCL_STATIC);
                return -1;
            }
        }
    } else {
        // do a double-to-int conversion first
        for (int i = 0; i < num; i++) {
            if (Tcl_GetDoubleFromObj(interp, data[i], &d) != TCL_OK) {
                Tcl_SetResult(interp, (char *) "Cgmap: error parsing vector element as integer", TCL_STATIC);
                return -1;
            }
            vec[i] = int (d);
        }
    }
    return num;
}

static int obj_Cgmap(ClientData /*UNUSED*/, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[])
{

    Tcl_Obj *atomselect = NULL;
    Tcl_Obj *object = NULL;
    Tcl_Obj *bytes = NULL;
    Tcl_Obj *bytes_append = NULL;
    Tcl_Obj *sel = NULL;
    float *coords = NULL;
    float *coords_append = NULL;
    const char *blockid_field = "user";
    const char *order_field = "user2";
    const char *weight_field= "user3";

    int nframes, natoms, ncoords, result, length;
    int first, last, stride;
    int molid, append_molid;

    natoms = ncoords = result = 0;
    molid = append_molid = 0;
    first = last = 0;
    stride = 1;
    nframes = 1;

    std::vector<float> weight;
    std::vector<int> bead;
    std::vector<int> index;

    // Parse Arguments
    int n = 1;
    while (n < argc) {
        const char *cmd = Tcl_GetString(objv[n]);
        if (!strncmp(cmd, "-molid", 7)) {
            if (Tcl_GetIntFromObj(interp,objv[n+1], &molid) != TCL_OK) {return TCL_ERROR;}
            n += 2;

        } else if (!strncmp(cmd, "-append", 8)) {
            if (Tcl_GetIntFromObj(interp,objv[n+1], &append_molid) != TCL_OK) {return TCL_ERROR;}
            n += 2;

        } else if (!strncmp(cmd, "-sel", 5)) {
            sel = objv[n+1];
            n += 2;

        } else if (!strncmp(cmd, "-first", 5)) {
            if (Tcl_GetIntFromObj(interp,objv[n+1], &first) != TCL_OK) {return TCL_ERROR;}
            n += 2;

        } else if (!strncmp(cmd, "-last", 4)) {
            if (Tcl_GetIntFromObj(interp,objv[n+1], &last) != TCL_OK) {return TCL_ERROR;}
            n += 2;

        } else if (!strncmp(cmd, "-stride", 6)) {
            if (Tcl_GetIntFromObj(interp,objv[n+1], &stride) != TCL_OK) {return TCL_ERROR;}
            n += 2;

        } else if (!strncmp(cmd, "-weight", 7)) {
            weight_field = Tcl_GetString(objv[n+1]);
            n += 2;

        } else if (!strncmp(cmd, "-blockid", 7)) {
            blockid_field = Tcl_GetString(objv[n+1]);
            n += 2;

        } else if (!strncmp(cmd, "-order", 6)) {
            order_field = Tcl_GetString(objv[n+1]);
            n += 2;

        } else {
            Tcl_WrongNumArgs(interp,1,objv, (char *)"molid");
            return TCL_ERROR;
        }
    }

    // Create an internal selection that we can manipulate if none was defined
    // Note that a passed selection overides the passed molid
    if (!sel) {
        Tcl_Obj *script = Tcl_ObjPrintf("atomselect %i all", molid);
        if (Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Cgmap: error calling atomselect", TCL_STATIC);
            return TCL_ERROR;
        }
        atomselect = Tcl_GetObjResult(interp);
        Tcl_IncrRefCount(atomselect);

    } else {
        // Create a internal selection that is a COPY of the passed selection
        atomselect = Tcl_DuplicateObj(sel);
        Tcl_IncrRefCount(atomselect);

        // Get the molid
        Tcl_Obj *script = Tcl_DuplicateObj(sel);
        Tcl_AppendToObj(script, " molid", -1);
        if(Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Cgmap: error calling atomselect", TCL_STATIC);
            return TCL_ERROR;
        }
        Tcl_Obj *molid_result =  Tcl_GetObjResult(interp);
        if (Tcl_GetIntFromObj(interp, molid_result, &molid) != TCL_OK) {return TCL_ERROR;}
    }

    // Get the number of frames
    Tcl_Obj *script = Tcl_ObjPrintf("molinfo %i get numframes", molid);
    if (Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error calling molinfo for nframes", TCL_STATIC);
        return TCL_ERROR;
    }
    object = Tcl_GetObjResult(interp);
    if (Tcl_GetIntFromObj(interp, object, &nframes) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing number of frames", TCL_STATIC);
        return TCL_ERROR;
    }

    if ( first < 0 || first >= nframes ) {
        Tcl_SetResult(interp, (char *) "Cgmap: illegal value of first_frame", TCL_STATIC);
        return TCL_ERROR;
    }
    if ( last == -1 || last > nframes || last < first ) last = nframes;

    // Get the number of atoms from selection
    script = Tcl_DuplicateObj(atomselect);
    Tcl_AppendToObj(script, " num", -1);
    if (Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error calling atomselect", TCL_STATIC);
        return TCL_ERROR;
    }
    object = Tcl_GetObjResult(interp);
    if (Tcl_GetIntFromObj(interp, object, &natoms) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing number of atoms", TCL_STATIC);
        return TCL_ERROR;
    }

    // Make sure we actually have some atoms
    if (natoms == 0) {
        Tcl_SetResult(interp, (char *) "Cgmap: Selection or molecule contains no atoms", TCL_STATIC);
        return TCL_ERROR;
    }

    // Get the weights (mass)
    script = Tcl_DuplicateObj(atomselect);
    Tcl_AppendPrintfToObj (script, " get %s", weight_field);
    if (Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error calling atomselect for weights", TCL_STATIC);
        return TCL_ERROR;
    }
    ncoords = parse_vector(Tcl_GetObjResult(interp), weight, interp);
    if (ncoords == -1) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing atomselect result", TCL_STATIC);
        return TCL_ERROR;
    }

    // Get the bead IDs
    script = Tcl_DuplicateObj(atomselect);
    Tcl_AppendPrintfToObj (script, " get %s", blockid_field);
    if (Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error calling atomselect for blocks", TCL_STATIC);
        return TCL_ERROR;
    }
    ncoords = parse_ivector(Tcl_GetObjResult(interp), bead, interp, true);
    if (ncoords == -1) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing atomselect result", TCL_STATIC);
        return TCL_ERROR;
    }

    // Get the atom IDs, we use these as a map when accessing the coordinate array
    // user2 is set via ::CGit::setBeadID
    script = Tcl_DuplicateObj(atomselect);
    Tcl_AppendPrintfToObj (script, " get %s", order_field);
    if (Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error calling atomselect for order", TCL_STATIC);
        return TCL_ERROR;
    }
    ncoords = parse_ivector(Tcl_GetObjResult(interp), index, interp, true);
    if (ncoords == -1) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing atomselect result", TCL_STATIC);
        return TCL_ERROR;
    }

    // Get current frame of the target mol
    script = Tcl_ObjPrintf("molinfo %d get frame", append_molid);
    if (Tcl_EvalObjEx(interp, script,  TCL_EVAL_DIRECT) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error getting append mol's current frame", TCL_STATIC);
        return TCL_ERROR;
    }
    int append_frame = 0;
    object = Tcl_GetObjResult(interp);
    if (Tcl_GetIntFromObj(interp, object, &append_frame) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing append mol's current frame", TCL_STATIC);
        return TCL_ERROR;
    }

    //Get number of atoms in target (append) mol
    script = Tcl_ObjPrintf("molinfo %i get numatoms", append_molid);
    if (Tcl_EvalObjEx(interp, script,  TCL_EVAL_DIRECT) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error getting append mol's number of atoms", TCL_STATIC);
        return TCL_ERROR;
    }
    int append_natoms = 0;
    object = Tcl_GetObjResult(interp);
    if (Tcl_GetIntFromObj(interp, object, &append_natoms) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "Cgmap: error parsing append mol's number of atoms", TCL_STATIC);
        return TCL_ERROR;
    }

    int print = ((last - first) / 10);
    if (print < 10) print = 10;
    if (print > 100) print = 100;

    //Loop over frames, calculate COMS, set coordinates in target mol
    for (int frame = first;
         frame <= last && frame < nframes;
         frame += stride) {

        if (frame % print == 0) {
            //Tcl_Obj *msg = Tcl_ObjPrintf ("puts \"Mapping frame %i\"", frame);
            Tcl_Obj *msg = Tcl_ObjPrintf ("vmdcon -info \"CGit> Mapping frame %i\"", frame);
            result = Tcl_EvalObjEx(interp, msg, TCL_EVAL_DIRECT);
            if (result != TCL_OK) { return TCL_ERROR; }
        }

        //Update the frames
        Tcl_Obj *mol_chgframe = Tcl_ObjPrintf ("molinfo top set frame %i", frame);
        if (Tcl_EvalObjEx(interp, mol_chgframe, TCL_EVAL_DIRECT) != TCL_OK)
            return TCL_ERROR;

        // Get the coordinates of the molecules in the reference mol
        Tcl_Obj *get_ts = Tcl_ObjPrintf("gettimestep %d %i", molid, frame);
        if (Tcl_EvalObjEx(interp, get_ts,  TCL_EVAL_DIRECT) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Cgmap: error getting coordinates", TCL_STATIC);
            return TCL_ERROR;
        }

        bytes = Tcl_GetObjResult(interp);
        Tcl_IncrRefCount(bytes);
        Tcl_InvalidateStringRep (bytes);
        coords = reinterpret_cast<float *> (Tcl_GetByteArrayFromObj(bytes, &length));

        /** Create a new frame for append_mol **/
        Tcl_ObjPrintf("animate dup %i", append_molid);
        if (Tcl_EvalObjEx(interp, get_ts,  TCL_EVAL_DIRECT) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Cgmap: error adding frame to append mol", TCL_STATIC);
            return TCL_ERROR;
        }
        append_frame++;

        Tcl_Obj *setframe = Tcl_ObjPrintf("molinfo %i set frame %i; display update", molid, frame);
        if (Tcl_EvalObjEx(interp, setframe,  TCL_EVAL_DIRECT) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Cgmap: error updating source frame", TCL_STATIC);
            return TCL_ERROR;
        }

        // Copy PBC conditions
        Tcl_Obj *setpbc = Tcl_ObjPrintf("molinfo %i set {a b c} [molinfo %i get {a b c}]", append_molid, molid);
        if (Tcl_EvalObjEx(interp, setpbc,  TCL_EVAL_DIRECT) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Cgmap: error updating PBC", TCL_STATIC);
            return TCL_ERROR;
        }

        // Get the coordinates of the molecules in the target (append) mol
        get_ts = Tcl_ObjPrintf("gettimestep %d %i", append_molid, append_frame);
        if (Tcl_EvalObjEx(interp, get_ts,  TCL_EVAL_DIRECT) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "Cgmap: error getting coordinates", TCL_STATIC);
            return TCL_ERROR;
        }

        bytes_append = Tcl_GetObjResult(interp);
        Tcl_IncrRefCount(bytes_append);
        Tcl_InvalidateStringRep(bytes_append);
        coords_append = reinterpret_cast<float *> (Tcl_GetByteArrayFromObj(bytes_append, &length));


        //loop over coordinates and beads, calculate COMs
        int current_bead, current_atom;
        current_bead = current_atom = 0;

        // Nested loop to work on each bead at a time
        float w,x,y,z;
        int j = 0;
        for (int start_atom = 0; start_atom < natoms; ) {
            current_bead = bead[start_atom];

            w = x = y = z = 0;
            // Calculate COM for each bead
            for ( current_atom = start_atom;
                  current_atom < natoms && bead[current_atom] == current_bead;
                  current_atom++) {

                //Lookup the atom index from the selection
                unsigned int idx = index[current_atom];
                float tw = weight[current_atom];

                w += tw;
                x += tw * coords[3*idx];
                y += tw * coords[3*idx+1];
                z += tw * coords[3*idx+2];
            }

            if (w == 0) {
                Tcl_SetResult(interp, (char *) "Cgmap: Bad weight can't total zero", TCL_STATIC);
                return TCL_ERROR;
            }

            // Insert calculated COMS into append_mols coordinate array
            // Need to figure out some kind of bounds checking here...
            coords_append[3 * j    ] = x / w;
            coords_append[3 * j + 1] = y / w;
            coords_append[3 * j + 2] = z / w;

            start_atom = current_atom;
            j++;
        } // bead loop

        // call rawtimestep to set byte array for append_mol
        Tcl_Obj *set_ts[5];

        set_ts[0] = Tcl_NewStringObj("rawtimestep", -1);
        set_ts[1] = Tcl_ObjPrintf("%d",append_molid);
        set_ts[2] = bytes_append;
        set_ts[3] = Tcl_NewStringObj("-frame", -1);
        set_ts[4] = Tcl_NewIntObj(append_frame);

        if (Tcl_EvalObjv (interp, 5, set_ts, 0) != TCL_OK)
            return TCL_ERROR;

        //Cleanup
        Tcl_DecrRefCount(bytes);
        Tcl_DecrRefCount(bytes_append);

    } // Frame loop

    //Cleanup
    Tcl_DecrRefCount(atomselect);

    Tcl_SetResult(interp, (char *) "", TCL_STATIC);
    return TCL_OK;
}

extern "C" {

    /* register the plugin with the tcl interpreters */
#if defined(CGMAPTCLDLL_EXPORTS) && defined( _WIN32 )
#  undef TCL_STORAGE_CLASS
#  define TCL_STORAGE_CLASS DLLEXPORT

#define WIN32_LEAN_AND_MEAN /* Exclude rarely-used stuff from Windows headers */
#include <windows.h>

    BOOL APIENTRY DllMain( HANDLE hModule,
                           DWORD  ul_reason_for_call,
                           LPVOID lpReserved )
    {
        return TRUE;
    }

    EXTERN int Cgmap_Init(Tcl_Interp *interp)

#else

            int Cgmap_Init(Tcl_Interp *interp)

#endif
    {

#if defined(USE_TCL_STUBS)
        if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL)
            return TCL_ERROR;
        if (Tcl_PkgRequire(interp, "Tcl", TCL_VERSION, 0) == NULL)
            return TCL_ERROR;
#endif

        if (Tcl_PkgProvide(interp, PACKAGE_NAME, PACKAGE_VERSION) != TCL_OK)
            return TCL_ERROR;

        Tcl_CreateObjCommand(interp, "::CGit::cgmap", obj_Cgmap,
                             (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
        return TCL_OK;
    }
}
