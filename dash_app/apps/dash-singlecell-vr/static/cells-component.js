/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('cells', {
    schema: {
        count: {type: 'number'},
        radius: {type: 'number'},
        scale: {type: 'number'},
        colors: {type: 'array'},
        positions: {type: 'array'}
    },
    init: function() {
        const {count, radius, scale, colors, positions} = this.data;  
        this.count = count;
        this.instanceColors = new Float32Array(count * 3);
        this.instanceColorsBase = new Float32Array(this.instanceColors.length);

        const material = this.material = new THREE.MeshLambertMaterial({ flatShading: true });

        const colorParsChunk = [
            'attribute vec3 instanceColor;',
            'varying vec3 vInstanceColor;',
            '#include <common>'
        ].join( '\n' );

        const instanceColorChunk = [
            '#include <begin_vertex>',
            '\tvInstanceColor = instanceColor;'
        ].join( '\n' );

        const fragmentParsChunk = [
            'varying vec3 vInstanceColor;',
            '#include <common>'
        ].join( '\n' );

        const colorChunk = [
            'vec4 diffuseColor = vec4( diffuse * vInstanceColor, opacity );'
        ].join( '\n' );

        material.onBeforeCompile = function ( shader ) {

            shader.vertexShader = shader.vertexShader
                .replace( '#include <common>', colorParsChunk )
                .replace( '#include <begin_vertex>', instanceColorChunk );

            shader.fragmentShader = shader.fragmentShader
                .replace( '#include <common>', fragmentParsChunk )
                .replace( 'vec4 diffuseColor = vec4( diffuse, opacity );', colorChunk );

        };

        this.el.setAttribute("id", "cells");
    },
    setMatrix: function ( pos, scaler ) {

        const position = new THREE.Vector3();
        const rotation = new THREE.Euler();
        const quaternion = new THREE.Quaternion();
        const scale = new THREE.Vector3();

        return function ( matrix ) {

            position.x = pos[0] * scaler;
            position.y = pos[1] * scaler;
            position.z = pos[2] * scaler;

            rotation.x = 0;
            rotation.y = 0;
            rotation.z = 0;

            quaternion.setFromEuler( rotation );

            scale.x = scale.y = scale.z = 1;

            matrix.compose( position, quaternion, scale );

        };

    },
    update: function( oldData ) {
        if (oldData.radius !== this.data.radius) {
            this.el.object3D.remove( this.mesh );

            this.geometry = new THREE.SphereBufferGeometry(this.data.radius);

            this.matrix = new THREE.Matrix4();
            this.mesh = new THREE.InstancedMesh( this.geometry, this.material, this.data.count );

            for ( let i = 0; i < this.data.count; i ++ ) {

                this.setMatrix(this.data.positions[i], this.data.scale)( this.matrix );
                this.mesh.setMatrixAt( i, this.matrix );
    
            }
            
            this.el.object3D.add( this.mesh );
        }
        if (oldData.colors !== this.data.colors) {
            const newColors = this.data.colors;
            const color = new THREE.Color();

            for ( let i = 0; i < this.count; i ++ ) {

                color.set(newColors[i]);
                color.toArray(this.instanceColors, i * 3);

            }

            this.instanceColorsBase.set(this.instanceColors);
            this.geometry.setAttribute( 'instanceColor', new THREE.InstancedBufferAttribute( new Float32Array( this.instanceColors ), 3 ) );
            this.geometry.setAttribute( 'instanceColorBase', new THREE.BufferAttribute( new Float32Array( this.instanceColorsBase ), 3 ) );
        } 
    },
    remove: function () {

        const meshes = [];

        scene.traverse( function ( object ) {

            if ( object.isMesh ) meshes.push( object );

        } );

        for ( let i = 0; i < meshes.length; i ++ ) {

            const mesh = meshes[ i ];
            mesh.material.dispose();
            mesh.geometry.dispose();

            this.el.object3D.remove( mesh );

        }

    }
});

