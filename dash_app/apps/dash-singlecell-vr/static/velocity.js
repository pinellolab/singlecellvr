/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('velocity', {
    schema: {
        count: {type: 'number'},
        radius: {type: 'number'},
        scale: {type: 'number'},
        colors: {type: 'array'},
        positions: {type: 'array'},
        endPositions: {type: 'array'}
    },
    init: function() {
        this.dummy = new THREE.Object3D();
        this.startObject = new THREE.Object3D();
        this.endObject = new THREE.Object3D();
        this.instanceColors = new Float32Array(this.data.count * 3);
        this.instanceColorsBase = new Float32Array(this.instanceColors.length);
       
        this.rotations = [];
        for ( var i = 0; i < this.data.count; i ++ ) {
            var x = this.data.positions[i][0] * this.data.scale;
            var y = this.data.positions[i][1] * this.data.scale;
            var z = this.data.positions[i][2] * this.data.scale;

            var xEnd = x + this.data.endPositions[i][0] * this.data.scale;
            var yEnd = y + this.data.endPositions[i][1] * this.data.scale;
            var zEnd = z + this.data.endPositions[i][2] * this.data.scale;

            var t = new TWEEN.Tween({x: x, y: y, z: z}).to({x: xEnd, y: yEnd, z: zEnd}, 5000).start()
            t.repeat(Infinity)
            const rotation = this.getDirection({'x':x,'y':y,'z':z}, 
                                          {'x':xEnd,'y':yEnd,'z':zEnd});
            this.rotations.push([rotation.x, rotation.y, rotation.z]);
        }

        let mesh;
        let geometry;
        let material;
        const el = this.el;
        var self = this;
        const loader = new THREE.GLTFLoader();
        loader.load("/assets/arrow_optimized/arrow_optimized.glb", function ( model ) {
            model.scene.traverse(node => {
                if (node.geometry && node.material) {
                    geometry = node.geometry;
                }
            })
            geometry.scale( .03, .03, .03 );
            geometry.computeVertexNormals();

            material = new THREE.MeshBasicMaterial({ flatShading: true });
            var colorParsChunk = [
                'attribute vec3 instanceColor;',
                'varying vec3 vInstanceColor;',
                '#include <common>'
            ].join( '\n' );
    
            var instanceColorChunk = [
                '#include <begin_vertex>',
                '\tvInstanceColor = instanceColor;'
            ].join( '\n' );
    
            var fragmentParsChunk = [
                'varying vec3 vInstanceColor;',
                '#include <common>'
            ].join( '\n' );
    
            var colorChunk = [
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

            mesh = new THREE.InstancedMesh(geometry, material, self.data.count);
            mesh.rotation.set(Math.PI/2,0,0);
            mesh.updateWorldMatrix(true, true);
            mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
            geometry = mesh.geometry;
            geometry.applyMatrix(mesh.matrixWorld);
            mesh.rotation.set(0,0,0);

            const color = new THREE.Color();    
            for ( let i = 0; i < self.data.count; i ++ ) {
                color.set(self.data.colors[i]);
                color.toArray(self.instanceColors, i * 3);
            }
            self.instanceColorsBase.set(self.instanceColors);
            mesh.geometry.setAttribute( 'instanceColor', new THREE.InstancedBufferAttribute( new Float32Array( self.instanceColors ), 3 ) );
            mesh.geometry.setAttribute( 'instanceColorBase', new THREE.BufferAttribute( new Float32Array( self.instanceColorsBase ), 3 ) );
            
            el.object3D.add(mesh);
            self.mesh = mesh

            self.tweens = TWEEN.getAll();
            for ( let index = 0; index < self.tweens.length; index ++ ) {
                let t = self.tweens[index];
                t.onUpdate((t) => {
                    self.dummy.position.set(t.x, t.y, t.z);
                    self.dummy.rotation.x = self.rotations[index][0];
                    self.dummy.rotation.y = self.rotations[index][1];
                    self.dummy.rotation.z = self.rotations[index][2];
                    self.dummy.updateMatrix();
                    self.mesh.setMatrixAt(index , self.dummy.matrix );
                    self.mesh.instanceMatrix.needsUpdate = true;
                })
            }
        } );
        this.el.setAttribute("id", "velocity");
    },
    setMatrix: function (start) {
        if (this.mesh) {
            for ( let i = 0; i < this.data.count; i += 1 ) {
                let x = this.data.positions[i][0] * this.data.scale;
                let y = this.data.positions[i][1] * this.data.scale;
                let z = this.data.positions[i][2] * this.data.scale;

                var xEnd = x + this.data.endPositions[i][0] * this.data.scale;
                var yEnd = y + this.data.endPositions[i][1] * this.data.scale;
                var zEnd = z + this.data.endPositions[i][2] * this.data.scale;

                if (start) {
                    this.dummy.position.set(x, y, z);
                } else {
                    this.dummy.position.set(xEnd, yEnd, zEnd);
                }
                this.dummy.rotation.x = this.rotations[i][0];
                this.dummy.rotation.y = this.rotations[i][1];
                this.dummy.rotation.z = this.rotations[i][2];
                this.dummy.updateMatrix();
                this.mesh.setMatrixAt(i , this.dummy.matrix );
            }
            this.mesh.instanceMatrix.needsUpdate = true;
            this.mesh.frustumCulled = false;
        }
    },
    update: function(oldData) {
        if (this.mesh) {
            if (JSON.stringify(oldData.colors) !== JSON.stringify(this.data.colors)) {
                if (this.mesh.geometry) {
                    const newColors = this.data.colors;
                    const color = new THREE.Color();

                    for ( let i = 0; i < this.data.count; i ++ ) {
                        color.set(newColors[i]);
                        color.toArray(this.instanceColors, i * 3);
                    }
                    this.instanceColorsBase.set(this.instanceColors);
                    this.mesh.geometry.setAttribute( 'instanceColor', new THREE.InstancedBufferAttribute( new Float32Array( this.instanceColors ), 3 ) );
                    this.mesh.geometry.setAttribute( 'instanceColorBase', new THREE.BufferAttribute( new Float32Array( this.instanceColorsBase ), 3 ) );
                }
            } 
            if (JSON.stringify(oldData.endPositions) !== JSON.stringify(this.data.endPositions) || 
                JSON.stringify(oldData.positions) !== JSON.stringify(this.data.positions)) {
                TWEEN.removeAll();
                for ( var i = 0; i < this.data.count; i ++ ) {
                    var x = this.data.positions[i][0] * this.data.scale;
                    var y = this.data.positions[i][1] * this.data.scale;
                    var z = this.data.positions[i][2] * this.data.scale;
        
                    var xEnd = x + this.data.endPositions[i][0] * this.data.scale;
                    var yEnd = y + this.data.endPositions[i][1] * this.data.scale;
                    var zEnd = z + this.data.endPositions[i][2] * this.data.scale;
                    
                    const rotation = this.getDirection({'x':x,'y':y,'z':z}, 
                                          {'x':xEnd,'y':yEnd,'z':zEnd});
                    this.rotations[i] = [rotation.x, rotation.y, rotation.z];

                    var t = new TWEEN.Tween({x: x, y: y, z: z}).to({x: xEnd, y: yEnd, z: zEnd}, 5000).start()
                    t.repeat(Infinity)
                }
                this.tweens = TWEEN.getAll();
                for ( let index = 0; index < this.tweens.length; index ++ ) {
                    let t = this.tweens[index];
                    t.onUpdate((t) => {
                        this.dummy.position.set(t.x, t.y, t.z);
                        this.dummy.rotation.x = this.rotations[index][0];
                        this.dummy.rotation.y = this.rotations[index][1];
                        this.dummy.rotation.z = this.rotations[index][2];
                        this.dummy.updateMatrix();
                        this.mesh.setMatrixAt(index , this.dummy.matrix );
                        this.mesh.instanceMatrix.needsUpdate = true;
                    })
                }
            }
        }
        this.setMatrix()
    },
    tick: function(time) {
        if (this.mesh) {
            TWEEN.update(time)
            this.mesh.instanceMatrix.needsUpdate = true;
            this.mesh.frustumCulled = false;
        }
    },
    getDirection: function(start, end) {
        this.startObject.position.set(start.x, start.y, start.z);
        this.endObject.position.set(end.x, end.y, end.z);
        this.startObject.lookAt(this.endObject.position);
        return this.startObject.rotation;
    }
});

