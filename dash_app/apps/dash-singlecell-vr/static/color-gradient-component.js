/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('color-gradient', {
    schema: {
        colors: {type: 'array'},
        maxLabel: {type: 'number', default: 1},
        minLabel: {type: 'number', default: 0},
        medianLabel: {type: 'number', default: 0.5},
        height: {type: 'number'},
        width: {type: 'number'},
        verticalOffset: {type: 'number', default: 0},
    },

    init: function () {

        // Container
        this.height = this.data.height;
        this.width = this.data.width;
        this.el.setAttribute('geometry', `primitive: plane; height: ${this.height}; width: ${this.width};`);
        this.el.setAttribute('material', `shader: flat; side:front; transparent: true; opacity: 0; alphaTest: 0.5; color: white;`);  

        // Scale Bar
        this.verticalOffset = -this.data.verticalOffset;
        this.colors = this.data.colors;
        const blockHeight = 5 / this.colors.length;
        const scaleBarX = this.width / 2 + .1;
        const scaleBarYBottom = this.verticalOffset - (blockHeight / 2);
        const scaleBarYTop = this.verticalOffset + (blockHeight * (this.colors.length - 1)) + (blockHeight / 2);
        const scaleBarZ = 0;
        const scaleBar = `<a-entity meshline="lineWidth: 5; path: ${scaleBarX} ${scaleBarYBottom} ${scaleBarZ}, ${scaleBarX} ${scaleBarYTop} ${scaleBarZ}; color: black"></a-entity>`;
        this.el.append(Utils.htmlToElement(scaleBar));
        const maxLabelText = `<a-text width="7" value="${this.data.maxLabel}" position="${scaleBarX + .1} ${scaleBarYTop} ${scaleBarZ}" color="black"></a-text>`;
        const minLabelText = `<a-text width="7" value="${this.data.minLabel}" position="${scaleBarX + .1} ${scaleBarYBottom} ${scaleBarZ}" color="black"></a-text>`;
        const medianLabelText = `<a-text width="7" value="${this.data.medianLabel}" position="${scaleBarX + .1} ${(scaleBarYTop + scaleBarYBottom) / 2} ${scaleBarZ}" color="black"></a-text>`;
        this.el.append(Utils.htmlToElement(maxLabelText), Utils.htmlToElement(minLabelText), Utils.htmlToElement(medianLabelText));
        
        // Instanced Blocks
        this.instanceColors = new Float32Array(this.colors.length * 3);
        this.instanceColorsBase = new Float32Array(this.instanceColors.length);

        const geometry = new THREE.PlaneBufferGeometry( this.width, blockHeight );
        this.geometry = geometry;
        const material = new THREE.MeshBasicMaterial({ flatShading: true });
        const positions = this.getPositions(this.colors.length, this.verticalOffset, blockHeight);

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

        const matrix = new THREE.Matrix4();
        const mesh = new THREE.InstancedMesh( geometry, material, this.colors.length );

        for ( let i = 0; i < this.colors.length; i ++ ) {
            this.setMatrix( positions[i] )( matrix );
            mesh.setMatrixAt( i, matrix );
        }
       
        this.el.object3D.add(mesh);
    },
    getPositions(num, verticalOffset, blockHeight) {
        positions = [];
        for (let i = 0; i < num; i++) {
            positions.push([0, verticalOffset + (blockHeight * i), 0]);
        }
        return positions;
    },
    setMatrix: function ( pos ) {

        const position = new THREE.Vector3();
        const rotation = new THREE.Euler();
        const quaternion = new THREE.Quaternion();
        const scale = new THREE.Vector3();

        return function ( matrix ) {

            position.x = pos[0];
            position.y = pos[1];
            position.z = pos[2];

            rotation.x = 0;
            rotation.y = 0;
            rotation.z = 0;

            quaternion.setFromEuler( rotation );

            scale.x = scale.y = scale.z = 1;

            matrix.compose( position, quaternion, scale );

        };

    },
    update: function( oldData ) {
        const newColors = this.data.colors;
        const color = new THREE.Color();

        for ( let i = 0; i < this.colors.length; i ++ ) {

            color.set(newColors[i]);
            color.toArray(this.instanceColors, i * 3);

        }

        this.instanceColorsBase.set(this.instanceColors);
        this.geometry.setAttribute( 'instanceColor', new THREE.InstancedBufferAttribute( new Float32Array( this.instanceColors ), 3 ) );
        this.geometry.setAttribute( 'instanceColorBase', new THREE.BufferAttribute( new Float32Array( this.instanceColorsBase ), 3 ) );
    },
});