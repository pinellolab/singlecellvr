/*global THREE, AFRAME, Utils, TWEEN*/

AFRAME.registerComponent('cells', {
    schema: {
        count: {type: 'number'},
        radius: {type: 'number'},
        scale: {type: 'number'},
        colors: {type: 'array'},
        positions: {type: 'array'},
        endPositions: {type: 'array'},
        textureSrc:{
            type:'string',
            default:"/assets/sphere.png"
        },
    },
    init: function() {
        const geometry = this.geometry = new THREE.BufferGeometry();

        this.vertices = [];
        for ( var i = 0; i < this.data.count; i ++ ) {
            var x = this.data.positions[i][0] * this.data.scale;
            var y = this.data.positions[i][1] * this.data.scale;
            var z = this.data.positions[i][2] * this.data.scale;

            var xEnd = x + this.data.endPositions[i][0] * this.data.scale;
            var yEnd = y + this.data.endPositions[i][1] * this.data.scale;
            var zEnd = z + this.data.endPositions[i][2] * this.data.scale;

            var t = new TWEEN.Tween({x: x, y: y, z: z}).to({x: xEnd, y: yEnd, z: zEnd}, 1000).start()
            t.repeat(Infinity)
            this.vertices.push( x, y, z );
        }
                
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
        geometry.setAttribute( 'position', new THREE.Float32BufferAttribute( this.vertices, 3 ) );
        geometry.computeBoundingSphere();
        this.boundingSphere = geometry.boundingSphere;
        this.points = new THREE.Points( geometry, material );
        this.el.object3D.add( this.points );
        this.el.setAttribute("id", "cells");
        
        this.tweens = TWEEN.getAll();
        for ( let index = 0; index < this.tweens.length; index ++ ) {
            let t = this.tweens[index];
            t.onUpdate((t) => {
                this.geometry.attributes.position.array[index*3] = t.x;
                this.geometry.attributes.position.array[index*3 + 1] = t.y;
                this.geometry.attributes.position.array[index*3 + 2] = t.z;
            })
        }
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

    },
    tick: function(time, timedelta) {
        TWEEN.update(time)
        // for ( var i = 0; i < this.tweens.length; i ++ ) {
        //     let x = this.tweens[i]._object.x;
        //     let y = this.tweens[i]._object.y;
        //     let z = this.tweens[i]._object.z;
            
        //     this.geometry.attributes.position.array[i*3] = x;
        //     this.geometry.attributes.position.array[i*3 + 1] = y;
        //     this.geometry.attributes.position.array[i*3 + 2] = z;
        // }
        this.geometry.attributes.position.needsUpdate = true;
    }
});