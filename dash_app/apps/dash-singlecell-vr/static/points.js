/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('cells', {
    schema: {
        count: {type: 'number'},
        radius: {type: 'number'},
        scale: {type: 'number'},
        colors: {type: 'array'},
        positions: {type: 'array'},
        textureSrc:{
            type:'string',
            default:"/assets/sphere.png"
        },
    },
    init: function() {

        const vertices = [];
        for ( var i = 0; i < this.data.count; i ++ ) {
            var x = this.data.positions[i][0] * this.data.scale;
            var y = this.data.positions[i][1] * this.data.scale;
            var z = this.data.positions[i][2] * this.data.scale;

            vertices.push( x, y, z );
        }

        const geometry = this.geometry = new THREE.BufferGeometry();
                
        const material = this.material = new THREE.PointsMaterial( { 
            size: this.data.radius,
            map: new THREE.TextureLoader().load(this.data.textureSrc),
            vertexColors: THREE.VertexColors,
            alphaTest: 0.5,
            transparent: false,
            depthTest: true,
            sizeAttenuation: true,
        });

        this.instanceColorsBase = new Float32Array(this.data.count * 3);
        geometry.setAttribute( 'color', this.instanceColorsBase )
        geometry.setAttribute( 'position', new THREE.Float32BufferAttribute( vertices, 3 ) );
        geometry.computeBoundingSphere();
        this.boundingSphere = geometry.boundingSphere;
        this.points = new THREE.Points( geometry, material );
        this.el.object3D.add( this.points );
        this.el.setAttribute("id", "cells");
    },
    update: function( oldData ) {
        const newColors = this.data.colors;
        const color = new THREE.Color();

        for ( let i = 0; i < this.data.count; i ++ ) {

            color.set(newColors[i]);
            color.toArray(this.instanceColorsBase, i * 3);

        }

        this.geometry.setAttribute( 'color', new THREE.BufferAttribute( new Float32Array( this.instanceColorsBase ), 3 ) );

        if (oldData.radius !== this.data.radius) {
            this.material.size = this.data.radius;
        }
    }
});