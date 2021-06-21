
#include "shader.hpp"

//--------------------------------------------------------------------------------------

oglu::shader::shader(const char vertex[], const char fragment[])
{
	/* 
		This constructor uses shaders stored in a header (*.hpp) file.
		Therefore compilation at runtime is not required.
		Usage.
		1) Include header file:
			#include "shadervars.hpp"
		2) Call constructor:
			oglu::shader cube_shader(cubeVertexShaderSource, cubeFragmentShaderSource);
	*/

	// Create the shaders
	unsigned int VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	unsigned int FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
		
	std::string FragmentShaderCode(fragment);
	std::string VertexShaderCode(vertex);
    
	char const * VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
	glCompileShader(VertexShaderID);

	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
	glCompileShader(FragmentShaderID);

	// Check compilation
	int status;
	int InfoLogLength;

	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &status);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	
	if ( InfoLogLength > 0 ){
		std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		printf("%s\n", &VertexShaderErrorMessage[0]);
	}

	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &status);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	
	if ( InfoLogLength > 0 ){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		printf("%s\n", &FragmentShaderErrorMessage[0]);
	}

	programID = glCreateProgram();
	
	glAttachShader(programID, VertexShaderID);
	glAttachShader(programID, FragmentShaderID);
	glLinkProgram(programID);

	// Check the program
	glGetProgramiv(programID, GL_LINK_STATUS, &status);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	
	if ( InfoLogLength > 0 ){
		std::vector<char> ProgramErrorMessage(InfoLogLength+1);
		glGetProgramInfoLog(programID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		printf("%s\n", &ProgramErrorMessage[0]);
	}
	
	glDetachShader(programID, VertexShaderID);
	glDetachShader(programID, FragmentShaderID);
	
	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);
}

//--------------------------------------------------------------------------------------

oglu::shader::shader(std::vector<const char *> &file_path)
{

	/*  
		This constructor reads shaders from disc during runtime.
		Usage.
		1) Prepare file paths:
			std::vector<const char *> shaders_file_path;
			shaders_file_path.push_back("shaders/CubeVertShader.vshdr");
			shaders_file_path.push_back("shaders/CubeFragmShader.fshdr");
		2) Call constructor:
			oglu::shader cube_shader(shaders_file_path);
	*/

	const char * vertex_file_path = file_path[0];
	const char * fragment_file_path = file_path[1];

	// Create the shaders
	unsigned int VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	unsigned int FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	std::ifstream VertexShaderStream;
	std::ifstream FragmentShaderStream;
	std::string FragmentShaderCode;
	std::string VertexShaderCode;

    // ensure ifstream objects can throw exceptions:
    VertexShaderStream.exceptions (std::ifstream::failbit | std::ifstream::badbit);
    FragmentShaderStream.exceptions (std::ifstream::failbit | std::ifstream::badbit);
    
	try 
    {
		// Read the Vertex Shader code from the file
		
		VertexShaderStream.open(vertex_file_path, std::ios::in);
		std::stringstream sstr1;
		sstr1 << VertexShaderStream.rdbuf();
		VertexShaderCode = sstr1.str();
		VertexShaderStream.close();

		// Read the Fragment Shader code from the file
		
		FragmentShaderStream.open(fragment_file_path, std::ios::in);
		std::stringstream sstr2;
		sstr2 << FragmentShaderStream.rdbuf();
		FragmentShaderCode = sstr2.str();
		FragmentShaderStream.close();
	}
    catch(std::ifstream::failure e)
    {
        std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
    }

	// Compile Vertex Shader
	printf("Compiling shader : %s\n", vertex_file_path);

	char const * VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
	glCompileShader(VertexShaderID);

	// Compile Fragment Shader
	printf("Compiling shader : %s\n", fragment_file_path);

	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
	glCompileShader(FragmentShaderID);

	// Check compilation
	int status;
	int InfoLogLength;

	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &status);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	
	if ( InfoLogLength > 0 ){
		std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		printf("%s\n", &VertexShaderErrorMessage[0]);
	}

	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &status);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	
	if ( InfoLogLength > 0 ){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		printf("%s\n", &FragmentShaderErrorMessage[0]);
	}

	// Link the program
	printf("Linking program\n");

	programID = glCreateProgram();
	
	glAttachShader(programID, VertexShaderID);
	glAttachShader(programID, FragmentShaderID);
	glLinkProgram(programID);

	// Check the program
	glGetProgramiv(programID, GL_LINK_STATUS, &status);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	
	if ( InfoLogLength > 0 ){
		std::vector<char> ProgramErrorMessage(InfoLogLength+1);
		glGetProgramInfoLog(programID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		printf("%s\n", &ProgramErrorMessage[0]);
	}
	
	glDetachShader(programID, VertexShaderID);
	glDetachShader(programID, FragmentShaderID);
	
	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);

}

//--------------------------------------------------------------------------------------

void oglu::shader::apply()
{
	glUseProgram(programID);
} 

//--------------------------------------------------------------------------------------

void oglu::shader::setUniformBool(const std::string &name, bool value) const
{         
    glUniform1i(glGetUniformLocation(programID, name.c_str()), (int)value); 
}

//--------------------------------------------------------------------------------------

void oglu::shader::setUniformInt(const std::string &name, int value) const
{ 
    glUniform1i(glGetUniformLocation(programID, name.c_str()), value); 
}

//--------------------------------------------------------------------------------------

void oglu::shader::setUniformFloat(const std::string &name, float value) const
{ 
    glUniform1f(glGetUniformLocation(programID, name.c_str()), value); 
}

//--------------------------------------------------------------------------------------

void oglu::shader::setUniformMatrFloat(const std::string &name, glm::mat4 value) const
{ 
    glUniformMatrix4fv(glGetUniformLocation(programID, name.c_str()), 1, GL_FALSE, glm::value_ptr(value)); 
}

//--------------------------------------------------------------------------------------

void oglu::shader::setUniformVecFloat(const std::string &name, glm::vec3 value) const
{ 
    glUniform3fv(glGetUniformLocation(programID, name.c_str()), 1, glm::value_ptr(value)); 
}

//--------------------------------------------------------------------------------------

