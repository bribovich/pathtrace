#include "scene.h"
#include "intersect.h"
#include "montecarlo.h"
#include "animation.h"

#include <thread>
using std::thread;

// modify the following line to disable/enable parallel execution of the pathtracer
bool parallel_pathtrace = true;

image3f pathtrace(Scene* scene, bool multithread);
void pathtrace(Scene* scene, image3f* image, RngImage* rngs, int offset_row, int skip_row, bool verbose);



// lookup texture value
vec3f lookup_scaled_texture(vec3f value, image3f* texture, vec2f uv, bool tile = false) {
    if(not texture) return value;
    
    float u=uv.x;
    float v=uv.y;

    int i= (int) (u*(texture->width()));
    int j=(int) (v*(texture->height()));

    int i_prime=i+1;
    int j_prime=j+1;

    float s = u*texture->width()-i;
    float t= v*texture->height()-j;

    if(tile){
        //u = ((int)(uv.x*(texture->width()-1)) % texture->width())/(float)texture->width();
        //v = (int)(uv.y*(texture->height()-1)) % texture->height()/(float)texture->height();
        i=i%(texture->width());
        if(i<0){
            i+=(texture->width());
        }
        i_prime=i_prime%(texture->width());
        if(i_prime<0){
            i_prime+=(texture->width());
        }

        j=j%(texture->height());
        if(j<0){
            j+=(texture->height());
        }

        j_prime=j_prime%(texture->height());
        if(j_prime<0){
            j_prime+=(texture->height());
        }



    }

    // for now, simply clamp texture coords
    else{
        i = clamp(i, 0, texture->width()-1);
        j = clamp(j, 0, texture->height()-1);

        i_prime = clamp(i_prime, 0, texture->width()-1);
        j_prime = clamp(j_prime, 0, texture->height()-1);

    }

    //return value *texture->at(i,j);
    return value * (texture->at(i, j)*(1.0f-s)*(1.0f-t)+texture->at(i_prime, j)*(s)*(1.0f-t)+texture->at(i, j_prime)*(1.0f-s)*(t)+texture->at(i_prime, j_prime)*(s)*(t));
}

// compute the brdf
vec3f eval_brdf(vec3f kd, vec3f ks, float n, vec3f v, vec3f l, vec3f norm, bool microfacet) {
    auto h = normalize(v+l);
    if (!microfacet) {

        return kd/pif + ks*(n+8)/(8*pif) * pow(max(0.0f,dot(norm,h)),n);
    } else {

        float d=(float(n+2.0f)/(2.0f*pif))*pow(max(0.0f,dot(norm,h)),n);
        vec3f f=ks+(one3f-ks)*pow(1.0f-dot(h,l),5);
        float g=min(min(1.0f,2.0*dot(h,norm)*dot(v,norm)/dot(v,h)), 2.0f*dot(h,norm)*dot(l,norm)/dot(l,h));
        return d*g*f/(4.0f*dot(l,norm)*dot(v,norm));

    }
}

// evaluate the environment map
vec3f eval_env(vec3f ke, image3f* ke_txt, vec3f dir) {
    //return ke; // <- placeholder
    float u = atan2(dir.x, dir.z)/(2.0f*pif);
    float v = 1.0f - acos(dir.y)/pif;
    return lookup_scaled_texture(ke, ke_txt, vec2f(u,v), true);

}

// compute the color corresponing to a ray by pathtrace
vec3f pathtrace_ray(Scene* scene, ray3f ray, Rng* rng, int depth) {
    // get scene intersection
    auto intersection = intersect(scene,ray);
    
    // if not hit, return background (looking up the texture by converting the ray direction to latlong around y)
    if((not (intersection.hit))|| (depth==scene->path_max_depth && depth >0)) {
        return eval_env(scene->background, scene->background_txt, ray.d);
    }

    
    // setup variables for shorter code
    auto pos = intersection.pos;
    auto norm = intersection.norm;
    auto v = -ray.d;
    
    // compute material values by looking up textures
    auto ke = lookup_scaled_texture(intersection.mat->ke, intersection.mat->ke_txt, intersection.texcoord, true);
    auto kd = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt, intersection.texcoord, true);
    auto ks = lookup_scaled_texture(intersection.mat->ks, intersection.mat->ks_txt, intersection.texcoord, true);
    auto n = intersection.mat->n;
    auto mf = intersection.mat->microfacet;
    
    // accumulate color starting with ambient
    auto c = scene->ambient * kd;
    
    // add emission if on the first bounce
    if(depth == 0 and dot(v,norm) > 0) c += ke;
    
    // foreach point light
    for(auto light : scene->lights) {
        // compute light response
        auto cl = light->intensity / (lengthSqr(light->frame.o - pos));
        // compute light direction
        auto l = normalize(light->frame.o - pos);
        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,l),0.0f) * eval_brdf(kd, ks, n, v, l, norm, mf);
        // multiply brdf and light
        auto shade = cl * brdfcos;
        // check for shadows and accumulate if needed
        if(shade == zero3f) continue;
        // if shadows are enabled
        if(scene->path_shadows) {
            // perform a shadow check and accumulate
            if(not intersect_shadow(scene,ray3f::make_segment(pos,light->frame.o))) c += shade;
        } else {
            // else just accumulate
            c += shade;
        }
    }
    
    // foreach surface
    for(Surface* s : scene->surfaces){

        // skip if no emission from surface
        if(s->mat->ke==zero3f){
            continue;
        }

        // todo: pick a point on the surface, grabbing normal, area, and texcoord
        // generate a 2d random number
        vec2f rand=rng->next_vec2f();
        float area;
        vec3f light_pos, light_norm;
        // check if quad
        if(s->isquad){

            // compute light position, normal, area

            light_pos=(s->frame.o+(rand.x-0.5)*2.0f*s->radius*s->frame.x+(rand.y-0.5)*2.0f*s->radius*s->frame.y);
            light_norm=normalize(s->frame.z);
            area=4.0f*(s->radius)*(s->radius);

            // set tex coords as random value got before

        }
        // else if sphere
        else{

            // compute light position, normal, area
            vec3f samp=normalize(sample_direction_spherical_uniform(rand));
            light_pos=s->frame.o+s->radius*samp;
            light_norm=normalize(samp);
            area=4.0f*pif*(s->radius)*(s->radius);
            // set tex coords as random value got before
        }
        intersection.texcoord=rand;
        // get light emission from material and texture
        vec3f kel=lookup_scaled_texture(s->mat->ke, s->mat->ke_txt, rand, true);
        // compute light direction
        vec3f light_direction=normalize(light_pos-pos);
        // compute light response (kel * area * cos_of_light / dist^2)
        vec3f light_response=kel*area*max(dot(light_norm, -light_direction),0.0f)/lengthSqr(light_pos-pos);
        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,light_direction),0.0f) * eval_brdf(kd, ks, n, v, light_direction, norm, mf);
        // multiply brdf and light
        auto response=brdfcos*light_response;
        // check for shadows and accumulate if needed
        // if shadows are enabled
        if(scene->path_shadows) {

            // perform a shadow check and accumulate
            if(not intersect_shadow(scene,ray3f::make_segment(pos, light_pos))) c+=response;
        } else {
            // else just accumulate
            return c+=response;

        }
    }
    
    // todo: sample the brdf for environment illumination if the environment is there
    // if scene->background is not zero3f
    if(scene->background!=zero3f){

        // pick direction and pdf

        pair<vec3f,float>result=sample_brdf(kd, ks, n, v, norm, rng->next_vec2f(),rng->next_float());
        vec3f d=result.first;
        float p=result.second;

        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,d),0.0f) * eval_brdf(kd, ks, n, v, d, norm, mf);
        // todo: accumulate response scaled by brdf*cos/pdf
        // if material response not zero3f
        if(brdfcos!=zero3f){
            // if shadows are enabled
            if(scene->path_shadows) {
                // perform a shadow check and accumulate
                if(not intersect_shadow(scene,ray3f(pos,d))) c+=brdfcos/p;
            } else {
                // else just accumulate
                c +=brdfcos/p;
            }
        }
    }
    // todo: sample the brdf for indirect illumination
    // if kd and ks are not zero3f and haven't reach max_depth
    if((depth < scene->path_max_depth) && (kd != zero3f || ks != zero3f)){
        // pick direction and pdf
        pair<vec3f,float>result=sample_brdf(kd, ks, n, v, norm, rng->next_vec2f(),rng->next_float());
        vec3f d=result.first;
        float p=result.second;
        // compute the material response (brdf*cos)
        vec3f brdfcos = max(dot(norm,d),0.0f) * eval_brdf(kd, ks, n, v, d, norm, mf);
        // accumulate recersively scaled by brdf*cos/pdf
        if(not (intersection.mat->ks == zero3f) ) {
        c= pathtrace_ray(scene, ray3f(pos, d),rng, depth+1)*brdfcos/p;
        }
        c += pathtrace_ray(scene, ray3f(pos, d),rng, depth+1)*brdfcos/p;
    }
    // if the material has reflections
    if(not (intersection.mat->kr == zero3f) ) {
        // create the reflection ray
        auto rr = ray3f(intersection.pos,reflect(ray.d,intersection.norm));
        // accumulate the reflected light (recursive call) scaled by the material reflection
        c += intersection.mat->kr * pathtrace_ray(scene,rr,rng,depth+1);
    }
    
    // return the accumulated color
    return c;
}


// runs the raytrace over all tests and saves the corresponding images
int main(int argc, char** argv) {
    auto args = parse_cmdline(argc, argv,
        { "04_pathtrace", "raytrace a scene",
            {  {"resolution",     "r", "image resolution", typeid(int),    true,  jsonvalue() } },
            {  {"scene_filename", "",  "scene filename",   typeid(string), false, jsonvalue("scene.json") },
               {"image_filename", "",  "image filename",   typeid(string), true,  jsonvalue("") } }
        });
    
    auto scene_filename = args.object_element("scene_filename").as_string();
    Scene* scene = nullptr;
    if(scene_filename.length() > 9 and scene_filename.substr(0,9) == "testscene") {
        int scene_type = atoi(scene_filename.substr(9).c_str());
        scene = create_test_scene(scene_type);
        scene_filename = scene_filename + ".json";
    } else {
        scene = load_json_scene(scene_filename);
    }
    error_if_not(scene, "scene is nullptr");
    
    auto image_filename = (args.object_element("image_filename").as_string() != "") ?
        args.object_element("image_filename").as_string() :
        scene_filename.substr(0,scene_filename.size()-5)+".png";
    
    if(not args.object_element("resolution").is_null()) {
        scene->image_height = args.object_element("resolution").as_int();
        scene->image_width = scene->camera->width * scene->image_height / scene->camera->height;
    }
    
    // NOTE: acceleration structure does not support animations
    message("reseting animation...\n");
    animate_reset(scene);
    
    message("accelerating...\n");
    accelerate(scene);
    
    message("rendering %s...\n", scene_filename.c_str());
    auto image = pathtrace(scene, parallel_pathtrace);
    
    message("saving %s...\n", image_filename.c_str());
    write_png(image_filename, image, true);
    
    delete scene;
    message("done\n");
}


/////////////////////////////////////////////////////////////////////
// Rendering Code


// pathtrace an image
void pathtrace(Scene* scene, image3f* image, RngImage* rngs, int offset_row, int skip_row, bool verbose) {
    if(verbose) message("\n  rendering started        ");
    // foreach pixel
    for(auto j = offset_row; j < scene->image_height; j += skip_row ) {
        if(verbose) message("\r  rendering %03d/%03d        ", j, scene->image_height);
        for(auto i = 0; i < scene->image_width; i ++) {
            // init accumulated color
            image->at(i,j) = zero3f;
            // grab proper random number generator
            auto rng = &rngs->at(i, j);
            // foreach sample
            for(auto jj : range(scene->image_samples)) {
                for(auto ii : range(scene->image_samples)) {
                    // compute ray-camera parameters (u,v) for the pixel and the sample
                    auto u = (i + (ii + rng->next_float())/scene->image_samples) /
                        scene->image_width;
                    auto v = (j + (jj + rng->next_float())/scene->image_samples) /
                        scene->image_height;
                    // compute camera ray
                    auto ray = transform_ray(scene->camera->frame,
                        ray3f(zero3f,normalize(vec3f((u-0.5f)*scene->camera->width,
                                                     (v-0.5f)*scene->camera->height,-1))));
                    // set pixel to the color raytraced with the ray
                    image->at(i,j) += pathtrace_ray(scene,ray,rng,0);
                }
            }
            // scale by the number of samples
            image->at(i,j) /= (scene->image_samples*scene->image_samples);
        }
    }
    if(verbose) message("\r  rendering done        \n");
    
}

// pathtrace an image with multithreading if necessary
image3f pathtrace(Scene* scene, bool multithread) {
    // allocate an image of the proper size
    auto image = image3f(scene->image_width, scene->image_height);
    
    // create a random number generator for each pixel
    auto rngs = RngImage(scene->image_width, scene->image_height);

    // if multitreaded
    if(multithread) {
        // get pointers
        auto image_ptr = &image;
        auto rngs_ptr = &rngs;
        // allocate threads and pathtrace in blocks
        auto threads = vector<thread>();
        auto nthreads = thread::hardware_concurrency();
        for(auto tid : range(nthreads)) threads.push_back(thread([=](){
            return pathtrace(scene,image_ptr,rngs_ptr,tid,nthreads,tid==0);}));
        for(auto& thread : threads) thread.join();
    } else {
        // pathtrace all rows
        pathtrace(scene, &image, &rngs, 0, 1, true);
    }
    
    // done
    return image;
}


