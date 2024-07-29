//convenience macros for error checking and attribute creation





#define MCHECKERROR(STAT, MSG) \
    if(STAT != MS::kSuccess){ \
		cerr << MSG << endl; \
		MGlobal::displayError(MSG); \
        return STAT; \
    }


#define MCHECKERROR_WARNING(STAT, MSG) \
	if(STAT != MS::kSuccess){ \
	cerr << MSG << endl; \
	MGlobal::displayWarning(MSG); \
	return STAT; \
	}


#define MCHECKERRORNORET(STAT, MSG) \
    if(STAT != MS::kSuccess){ \
		cerr << MSG << endl; \
		MGlobal::displayError(MSG); \
    }


#define MAKE_TYPED_ATTR(ATT_VARIABLE_NAME, TYPE, DEFAULT_VALUE) \
	MStatus ATT_VARIABLE_NAME##_stat; \
	MFnTypedAttribute ATT_VARIABLE_NAME##_fn; \
	ATT_VARIABLE_NAME = ATT_VARIABLE_NAME##_fn.create(#ATT_VARIABLE_NAME, #ATT_VARIABLE_NAME, TYPE, DEFAULT_VALUE); \
	ATT_VARIABLE_NAME##_fn.setHidden(true); \
	ATT_VARIABLE_NAME##_stat = addAttribute(ATT_VARIABLE_NAME); \
	MCHECKERROR(ATT_VARIABLE_NAME##_stat, "Error in adding new typed attribute");

#define MAKE_NUMERIC_ATTR(ATT_VARIABLE_NAME, TYPE, DEFAULT_VALUE, KEYABLE) \
	MStatus ATT_VARIABLE_NAME##_stat; \
	MFnNumericAttribute ATT_VARIABLE_NAME##_fn; \
	ATT_VARIABLE_NAME = ATT_VARIABLE_NAME##_fn.create(#ATT_VARIABLE_NAME, #ATT_VARIABLE_NAME, TYPE, DEFAULT_VALUE); \
	MCHECKERROR(ATT_VARIABLE_NAME##_stat, "numeric attr create error"); \
	ATT_VARIABLE_NAME##_fn.setKeyable(KEYABLE); \
	ATT_VARIABLE_NAME##_stat = addAttribute(ATT_VARIABLE_NAME); \
	MCHECKERROR(ATT_VARIABLE_NAME##_stat, "Error in adding new numeric attribute");

#define ADD_ATTRIBUTE(ATTR) \
	MStatus ATTR##_stat; \
	ATTR##_stat = addAttribute(ATTR); \
    MCHECKERROR(ATTR##_stat, "addAttribute: ATTR")

#define ATTRIBUTE_AFFECTS(IN, OUT) \
	MStatus IN##OUT##_stat; \
	IN##OUT##_stat = attributeAffects(IN, OUT); \
	MCHECKERROR(IN##OUT##_stat,"attributeAffects:" #IN "->" #OUT);



#define REGISTER_NODE(PLUGIN, NODE_TYPE, CLASS_NAME, CLASSIFICATION) \
	MStatus CLASS_NAME##_stat; \
	CLASS_NAME##_stat = PLUGIN.registerNode(CLASS_NAME::mTypeName , CLASS_NAME::mTypeId, CLASS_NAME::creator, CLASS_NAME::initialize, NODE_TYPE, CLASSIFICATION); \
	if(CLASS_NAME##_stat != MS::kSuccess) \
	{ \
		MString errorMsg("Failed to register user node " #CLASS_NAME); \
		cerr << errorMsg << endl; \
		MGlobal::displayError(errorMsg); \
		return CLASS_NAME##_stat; \
	} \




#define REGISTER_COMMAND(PLUGIN, CLASS_NAME) \
	MStatus CLASS_NAME##_stat; \
	CLASS_NAME##_stat = PLUGIN.registerCommand(CLASS_NAME::commandName() ,CLASS_NAME::creator, NULL); \
	if(CLASS_NAME##_stat != MS::kSuccess) \
	{ \
		MString errorMsg("Failed to register user command: " + CLASS_NAME::commandName()); \
		cerr << errorMsg << endl; \
		MGlobal::displayError(errorMsg); \
		return CLASS_NAME##_stat; \
	} \



#define REGISTER_COMMAND_WITH_SYNTAX(PLUGIN, CLASS_NAME) \
	MStatus CLASS_NAME##_stat; \
	CLASS_NAME##_stat = PLUGIN.registerCommand(CLASS_NAME::commandName() ,CLASS_NAME::creator, CLASS_NAME::syntax); \
	if(CLASS_NAME##_stat != MS::kSuccess) \
	{ \
		MString errorMsg("Failed to register user command: " + CLASS_NAME::commandName()); \
		cerr << errorMsg << endl; \
		MGlobal::displayError(errorMsg); \
		return CLASS_NAME##_stat; \
	} \


#define DEREGISTER_NODE(PLUGIN, CLASS_NAME) \
	MStatus CLASS_NAME##_stat; \
	CLASS_NAME##_stat = PLUGIN.deregisterNode(CLASS_NAME::mTypeId); \
	if(CLASS_NAME##_stat != MS::kSuccess) \
	{ \
		MString errorMsg("Failed to deregister user node: " #CLASS_NAME); \
		cerr << errorMsg << endl; \
		MGlobal::displayError(errorMsg); \
		return CLASS_NAME##_stat; \
	} \




#define DEREGISTER_COMMAND(PLUGIN, CLASS_NAME) \
	MStatus CLASS_NAME##_stat; \
	CLASS_NAME##_stat = PLUGIN.deregisterCommand(CLASS_NAME::commandName()); \
	if(CLASS_NAME##_stat != MS::kSuccess) \
	{ \
		MString errorMsg("Failed to deregister user command: " + CLASS_NAME::commandName()); \
		cerr << errorMsg << endl; \
		MGlobal::displayError(errorMsg); \
		return CLASS_NAME##_stat; \
	} \




