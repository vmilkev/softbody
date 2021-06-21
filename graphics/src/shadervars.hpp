#ifndef SHADERVARS_HPP
#define SHADERVARS_HPP

const char VertexShaderSource[] = R"glsl(
#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexCoord;
layout (location = 2) in vec3 aNormal;

uniform mat4 Transform;

uniform mat4 Model;
uniform mat4 View;
uniform mat4 Projection;
uniform mat4 camera;

out vec4 vertexColour;
out vec2 TexCoord;
out vec3 Normal;
out vec3 FragPos;

void main()
{
    gl_Position = Projection * View * Model * Transform * vec4(aPos, 1.0);
    FragPos = vec3(Model * vec4(aPos, 1.0));
    Normal = aNormal;  
    vertexColour = vec4(0.5, 0.0, 0.0, 1.0);
    TexCoord = aTexCoord;
}

)glsl";

const char colorVertexShaderSource[] = R"glsl(
#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

uniform mat4 Transform;

uniform mat4 Model;
uniform mat4 View;
uniform mat4 Projection;

out vec4 vertexColour;
out vec3 Normal;
out vec3 FragPos;

void main()
{
    gl_Position = Projection * View * Model * Transform * vec4(aPos, 1.0);
    FragPos = vec3(Model * vec4(aPos, 1.0));
    Normal = mat3(transpose(inverse(Model))) * aNormal;  
    vertexColour = vec4(0.5, 0.0, 0.0, 1.0);
}

)glsl";

const char FragmentShaderSource[] = R"glsl(
#version 330 core

in vec2 TexCoord;
in vec3 Normal;
in vec3 FragPos;

out vec4 FragColor;

uniform sampler2D texture_1;

uniform vec3 lightColor;
uniform vec3 lightPos;
uniform vec3 viewPos;

void main()
{
    // ambient
    float ambientStrength = 0.5;
    vec3 ambient = ambientStrength * lightColor * texture(texture_1, TexCoord).rgb;
  	
    // diffuse 
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor * texture(texture_1, TexCoord).rgb;

    // spectacular light calculations
    float specularStrength = 0.5;
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 64);
    vec3 specular = specularStrength * spec * lightColor * texture(texture_1, TexCoord).rgb; 

    vec3 result = (ambient + diffuse + specular);

    FragColor = vec4(result, 1.0);

	//FragColor = texture(texture_1, TexCoord);
}

)glsl";

const char colourFragmentShaderSource[] = R"glsl(
#version 330 core

in vec3 Normal;
in vec3 FragPos;

out vec4 FragColor;

uniform vec3 objectColor;
uniform vec3 lightColor;
uniform vec3 lightPos;
uniform vec3 viewPos;

void main()
{
    // ambient
    float ambientStrength = 0.3;
    vec3 ambient = ambientStrength * lightColor;
  	
    // diffuse 
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.1);
    vec3 diffuse = diff * lightColor;

    // spectacular light calculations
    float specularStrength = 0.5;
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 64);
    vec3 specular = specularStrength * spec * lightColor; 

    vec3 result = (ambient + diffuse + specular) * objectColor;

    FragColor = vec4(result, 1.0);
}

)glsl";

const char lightFragmentShaderSource[] = R"glsl(
#version 330 core

out vec4 FragColor;

uniform vec3 lightColor;

void main()
{
    FragColor = vec4(lightColor, 1.0f);
}

)glsl";

#endif